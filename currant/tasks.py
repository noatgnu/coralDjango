import os.path
from json import dumps, loads

import numpy as np
import pandas as pd
from alphastats import GenericLoader, DataSet
from asgiref.sync import async_to_sync
from channels.layers import get_channel_layer
from django_rq import job
from django.core.files.base import File as djangoFile, ContentFile
from currant.models import File, Operation
from currant.serializers import FileSerializer
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

from currant.utility import replace_special_with_dot, read_fasta, reg_positon_residue


@job('default')
def save_file(request):
    f = File(file_type="raw_file")
    channel_layer = get_channel_layer()
    async_to_sync(channel_layer.group_send)("1", {
        'type': 'event_message',
        'message': {
            'message': "Start uploaded",
            'senderName': "Server",
            'requestType': "fileUpload"
        }
    })
    extension = request.data['file'].name.split(".")[-1]
    f.file.save(f"{str(f.link_id)}.{extension}", djangoFile(request.data['file']))
    if extension == "csv":
        df = pd.read_csv(f.file.url)
    elif extension == "xlsx":
        df = pd.read_excel(f.file.url)
    elif extension == "txt" or extension == "tsv":
        df = pd.read_csv(f.file.url, sep="\t")
    df.replace("#VALUE!", np.NAN, inplace=True)
    f.columns = dumps(df.columns.tolist())
    f.save()
    file_json = FileSerializer(f, many=False, context={"request": request})
    async_to_sync(channel_layer.group_send)("1", {
        'type': 'event_message',
        'message': {
            'message': "Finished uploaded",
            'data': file_json.data,
            'senderName': "Server",
            'requestType': "fileUpload"
        }
    })

@job('default')
def r_qfeatures_protein(operation_id: int, session_id: str):
    ro.r("rm(list = ls(all.names = TRUE))")
    o = Operation.objects.get(id=operation_id)
    channel_layer = get_channel_layer()
    message_template = {
        'message': "Convert dataframe to R dataframe",
        'senderName': "Server",
        'requestType': "QFeatures (Protein) (R)",
        'operationId': o.id
    }
    request_form = loads(o.value)
    input_file_id = o.input_files.all()[0].id

    if request_form["operationType"] == "RQF-PEP":
        message_template["requestType"] = "QFeatures (Peptide) (R)"

    try:
        message_template["message"] = "Running operation"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        file = File.objects.get(id=input_file_id)
        extension = file.file.name.split(".")[-1]

        message_template["message"] = "Loading file"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        if extension == "csv":
            df = pd.read_csv(file.file.url)
        elif extension == "xlsx":
            df = pd.read_excel(file.file.url)
        elif extension == "txt" or extension == "tsv":
            df = pd.read_csv(file.file.url, sep="\t")
        df.replace("#VALUE!", np.NAN, inplace=True)

        df[request_form['sampleColumns']] = df[request_form['sampleColumns']].astype(float)
        # convert to R dataframe
        message_template["message"] = "Convert dataframe to R dataframe"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        with (ro.default_converter + pandas2ri.converter).context():
            r_from_pd_df = ro.conversion.get_conversion().py2rpy(df)
            qfeat = importr("QFeatures")
            limma = importr("limma")
        sampleColumnsIndex = [n+1 for n, c in enumerate(df.columns) if c in request_form['sampleColumns']]
        conditions = []
        samples = []
        for c in df.columns:
            if c in request_form['conditionMap']:
                cond = replace_special_with_dot(request_form['conditionMap'][c]["condition"])
                if cond[0].isdigit():
                    cond = "X"+cond
                conditions.append(cond)
                samples.append(replace_special_with_dot(request_form['conditionMap'][c]["replicate"]))
        sampleColumnsIndex = ro.IntVector(sampleColumnsIndex)
        ro.globalenv["conditions"] = ro.StrVector(conditions)
        ro.globalenv["samples"] = ro.StrVector(samples)
        message_template["message"] = "Creating QFeatures object"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })

        data = qfeat.readQFeatures(r_from_pd_df, ecol=sampleColumnsIndex, name='startingDF')
        ro.globalenv["data"] = data
        ro.r("""
        data$group <- conditions
        data$sample <- samples
        currentAssayName <- 'startingDF'
        """)
        message_template["message"] = "Filter NA"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        ro.r(f"""
        data <- selectRowData(data, c("{replace_special_with_dot(request_form["indexColumn"])}"))
        data <- zeroIsNA(data, i = seq_along(data))
        data <- filterNA(data, i = seq_along(data), pNA = {request_form["dataCompleteness"]})
        """)
        if request_form["imputation"]:
            message_template["message"] = "Perform imputation"
            async_to_sync(channel_layer.group_send)(session_id, {
                'type': 'job_message',
                'message': message_template
            })
            ro.r(f"""
            data <- impute(data, method = "{request_form["imputation"]}", i ="startingDF")
            currentAssayName <- "imputedAssay"
            """)
            d = ro.r("data")

        if request_form["log2"]:
            message_template["message"] = "Perform log2 transformation"
            async_to_sync(channel_layer.group_send)(session_id, {
                'type': 'job_message',
                'message': message_template
            })
            ro.r(f"""
            data <- addAssay(data, logTransform(data[[seq_along(data)[length(seq_along(data))]]]), name="log2")
            currentAssayName <- "log2"
            """)

        if request_form["normalization"]:
            message_template["message"] = "Perform normalization"
            async_to_sync(channel_layer.group_send)(session_id, {
                'type': 'job_message',
                'message': message_template
            })
            ro.r(f"""
            data <- addAssay(data, normalize(data[[seq_along(data)[length(seq_along(data))]]], method="{request_form["normalization"]}"), name="norm")
            currentAssayName <- "norm"
            """)

        if request_form["operationType"] == "RQF-PEP":
            message_template["message"] = "Aggregate features"
            async_to_sync(channel_layer.group_send)(session_id, {
                'type': 'job_message',
                'message': message_template
            })
            ro.r(f"""
            data <- aggregateFeatures(data, i = currentAssayName, fcol = {replace_special_with_dot(request_form["aggregateColumn"])}, name = "proteins")
            currentAssayName <- "proteins"
            """)

        conditionA = replace_special_with_dot(request_form["conditionA"])
        if conditionA[0].isdigit():
            conditionA = "X"+conditionA
        conditionB = replace_special_with_dot(request_form["conditionB"])
        if conditionB[0].isdigit():
            conditionB = "X"+conditionB
        message_template["message"] = "Create contrast matrix"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        ro.r(f"""
        design <- model.matrix(~0+data$group)
        colnames(design) <- gsub("data\\\\$group", "", colnames(design))
        fit <- lmFit(assay(data, currentAssayName), design)
    
        """)
        design = ro.r("design")
        conMat = limma.makeContrasts(contrast1=conditionA+"-"+conditionB, levels=design)
        ro.globalenv["contrast.matrix"] = conMat
        ro.r("""
        fit <- contrasts.fit(fit, contrast.matrix)
        fit <- eBayes(fit)
        """)
        message_template["message"] = "Perform differential analysis using limma"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        ro.r("""
        result <- topTable(fit, coef=1, adjust="BH", number=Inf, sort.by="none")
        result <- cbind(as.data.frame(rowData(data[[currentAssayName]])), result)
        """)

        result = ro.r("result")
        with (ro.default_converter + pandas2ri.converter).context():
            pd_from_r_df = ro.conversion.get_conversion().rpy2py(result)
        f = save_df(pd_from_r_df, "differential_analysis;table")
        o.output_files.set([f])
        o.job_finished = True
        o.save()
        message_template["message"] = "Completed operation"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
    except Exception as error:
        o.job_error = str(error)
        o.job_error_status = True
        o.save()
        message_template["message"] = "Error: " + str(error)
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })


@job('default')
def alphapeptstats_diff(operation_id: int, session_id: str):

    o = Operation.objects.get(id=operation_id)
    message_template = {
        'message': "Running operation",
        'senderName': "Server",
        'requestType': "Differential Analysis (AlphaPeptStats)",
        'operationId': o.id
    }
    channel_layer = get_channel_layer()
    try:
        form = loads(o.value)
        df = load_dataframe(channel_layer, o.input_files.all()[0], message_template, form, session_id)
        message_template["message"] = "Loaded input file"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        alpha_data = GenericLoader(file=df, intensity_column=form["sampleColumns"],
                                   index_column=form["indexColumn"])
        meta_df = []
        for i in form["conditionMap"]:
            meta_df.append({"sample": i, "condition": form["conditionMap"][i]["condition"]})
        meta_df = pd.DataFrame(meta_df)
        data_set = DataSet(
            loader=alpha_data,
            metadata_path=meta_df,
            sample_column="sample",
        )
        message_template["message"] = "Created dataset"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        data_set.preprocess(
            remove_contaminations=form["removeContaminants"],
            log2_transform=form["log2"],
            normalization=form["normalization"],
            imputation=form["imputation"],
            data_completeness=form["dataCompleteness"],
        )
        message_template["message"] = "Preprocessed dataset"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        result = data_set.diff_expression_analysis(form["conditionA"],
                                                   form["conditionB"],
                                                   "condition",
                                                   form["diffTest"])
        message_template["message"] = "Performed differential analysis"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        result = result.merge(df, on=form["indexColumn"], how="left")
        f = save_df(result, "differential_analysis;table")
        o.output_files.set([f])
        message_template["message"] = "Completed operation"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        o.job_finished = True
        o.save()
    except Exception as error:
        o.job_error = str(error)
        o.job_error_status = True
        o.save()
        message_template["message"] = "Error: " + str(error)
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })


@job('default')
def r_qfeatures_imputation(operation_id: int, session_id: str):
    ro.r("rm(list = ls(all.names = TRUE))")
    o = Operation.objects.get(id=operation_id)
    message_template = {
        'message': "Running operation",
        'senderName': "Server",
        'requestType': "QFeatures (Imputation) (R)",
        'operationId': o.id
    }
    channel_layer = get_channel_layer()
    request_form = loads(o.value)
    input_file_id = o.input_files.all()[0].id
    async_to_sync(channel_layer.group_send)(session_id, {
        'type': 'job_message',
        'message': message_template
    })
    try:
        df = load_dataframe(channel_layer, input_file_id, message_template, request_form, session_id)
        with (ro.default_converter + pandas2ri.converter).context():
            r_from_pd_df = ro.conversion.get_conversion().py2rpy(df)
            qfeat = importr("QFeatures")
            limma = importr("limma")
        sampleColumnsIndex = [n + 1 for n, c in enumerate(df.columns) if c in request_form['sampleColumns']]
        conditions = []
        samples = []
        for c in df.columns:
            if c in request_form['conditionMap']:
                cond = replace_special_with_dot(request_form['conditionMap'][c]["condition"])
                if cond[0].isdigit():
                    cond = "X" + cond
                conditions.append(cond)
                samples.append(replace_special_with_dot(request_form['conditionMap'][c]["replicate"]))
        sampleColumnsIndex = ro.IntVector(sampleColumnsIndex)
        ro.globalenv["conditions"] = ro.StrVector(conditions)
        ro.globalenv["samples"] = ro.StrVector(samples)
        message_template["message"] = "Creating QFeatures object"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })

        data = qfeat.readQFeatures(r_from_pd_df, ecol=sampleColumnsIndex, name='startingDF')
        ro.globalenv["data"] = data
        ro.r("""
        data$group <- conditions
        data$sample <- samples
        currentAssayName <- 'startingDF'
        """)
        message_template["message"] = "Filter NA"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        ro.r(f"""
        data <- selectRowData(data, c("{replace_special_with_dot(request_form["indexColumn"])}"))
        data <- zeroIsNA(data, i = seq_along(data))
        d <- filterNA(data, i = seq_along(data), pNA = {request_form["dataCompleteness"]})
        """)
        if request_form["imputation"]:
            message_template["message"] = "Perform imputation"
            async_to_sync(channel_layer.group_send)(session_id, {
                'type': 'job_message',
                'message': message_template
            })
            ro.r(f"""
            d <- impute(d, method = "{request_form["imputation"]}", i ="startingDF")
            currentAssayName <- "imputedAssay"
            """)
        ro.r("""
        result <- cbind(as.data.frame(rowData(d[["startingDF"]])), as.data.frame(assay(d, "startingDF")))
        """)

        result = ro.r("result")
        with (ro.default_converter + pandas2ri.converter).context():
            pd_from_r_df = ro.conversion.get_conversion().rpy2py(result)
        f = save_df(pd_from_r_df, "imputation;table")
        o.output_files.set([f])
        o.job_finished = True
        o.save()
        message_template["message"] = "Completed operation"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        ro.r("rm(list = ls(all.names = TRUE))")
    except Exception as error:
        o.job_error = str(error)
        o.job_error_status = True
        o.save()
        message_template["message"] = "Error: " + str(error)
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })


@job('default')
def r_qfeatures_normalization(operation_id: int, session_id: str):
    o = Operation.objects.get(id=operation_id)
    message_template = {
        'message': "Running operation",
        'senderName': "Server",
        'requestType': "QFeatures (Normalization) (R)",
        'operationId': o.id
    }
    channel_layer = get_channel_layer()
    request_form = loads(o.value)
    input_file_id = o.input_files.all()[0].id
    try:
        df = load_dataframe(channel_layer, input_file_id, message_template, request_form, session_id)
        with (ro.default_converter + pandas2ri.converter).context():
            r_from_pd_df = ro.conversion.get_conversion().py2rpy(df)
            qfeat = importr("QFeatures")
        sampleColumnsIndex = [n + 1 for n, c in enumerate(df.columns) if c in request_form['sampleColumns']]
        conditions = []
        samples = []
        for c in df.columns:
            if c in request_form['conditionMap']:
                cond = replace_special_with_dot(request_form['conditionMap'][c]["condition"])
                if cond[0].isdigit():
                    cond = "X" + cond
                conditions.append(cond)
                samples.append(replace_special_with_dot(request_form['conditionMap'][c]["replicate"]))
        sampleColumnsIndex = ro.IntVector(sampleColumnsIndex)
        ro.globalenv["conditions"] = ro.StrVector(conditions)
        ro.globalenv["samples"] = ro.StrVector(samples)
        message_template["message"] = "Creating QFeatures object"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })

        data = qfeat.readQFeatures(r_from_pd_df, ecol=sampleColumnsIndex, name='startingDF')
        ro.globalenv["data"] = data
        ro.r("""
                data$group <- conditions
                data$sample <- samples
                currentAssayName <- 'startingDF'
                """)
        message_template["message"] = "Perform normalization"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        ro.r(f"""
        data <- selectRowData(data, c("{replace_special_with_dot(request_form["indexColumn"])}"))
        """)
        if request_form["log2"]:
            message_template["message"] = "Perform log2 transformation"
            async_to_sync(channel_layer.group_send)(session_id, {
                'type': 'job_message',
                'message': message_template
            })
            ro.r(f"""
            data <- addAssay(data, logTransform(data[[currentAssayName]]), name="log2")
            currentAssayName <- "log2"
            """)
        result = ro.r(
            """cbind(as.data.frame(rowData(data[[currentAssayName]])), as.data.frame(assay(data, currentAssayName)))""")
        with (ro.default_converter + pandas2ri.converter).context():
            pd_from_r_df = ro.conversion.get_conversion().rpy2py(result)
        plotly_graph_data_before_normalization = {}
        for s in request_form["sampleColumns"]:
            if s in request_form["conditionMap"]:
                condition = request_form["conditionMap"][s]["condition"]
                if condition not in plotly_graph_data_before_normalization:
                    plotly_graph_data_before_normalization[condition] = {
                        "x": [],
                        "y": [],
                        "boxpoints": False,
                        "type": 'box',
                        "showlegend": False,
                        "name": condition
                    }
                value = list(pd_from_r_df[s].tolist())
                plotly_graph_data_before_normalization[condition]["x"] = plotly_graph_data_before_normalization[condition]["x"] + [s] * len(value)
                plotly_graph_data_before_normalization[condition]["y"] = plotly_graph_data_before_normalization[condition]["y"] + value

        ro.r(f"""
        data <- addAssay(data, normalize(data[[currentAssayName]], method="{request_form["normalization"]}"), name="norm")
        currentAssayName <- "norm"
        """)
        result = ro.r("""cbind(as.data.frame(rowData(data[[currentAssayName]])), as.data.frame(assay(data, currentAssayName)))""")

        with (ro.default_converter + pandas2ri.converter).context():
            pd_from_r_df = ro.conversion.get_conversion().rpy2py(result)
        plotly_graph_data_after_normalization = {}
        for s in request_form["sampleColumns"]:
            if s in request_form["conditionMap"]:
                condition = request_form["conditionMap"][s]["condition"]
                if condition not in plotly_graph_data_after_normalization:
                    plotly_graph_data_after_normalization[condition] = {
                        "x": [],
                        "y": [],
                        "boxpoints": False,
                        "type": 'box',
                        "showlegend": False,
                        "name": condition
                    }
                value = list(pd_from_r_df[s].tolist())
                plotly_graph_data_after_normalization[condition]["x"] = \
                plotly_graph_data_after_normalization[condition]["x"] + [s] * len(value)
                plotly_graph_data_after_normalization[condition]["y"] = \
                plotly_graph_data_after_normalization[condition]["y"] + value

        tick_val = []
        tick_text = []
        graph_box_before_normalization = []
        graph_box_after_normalization = []
        for c in plotly_graph_data_before_normalization:
            tick_val.append(plotly_graph_data_before_normalization[c]['x'][np.round(len(plotly_graph_data_before_normalization[c]['x'])/2).astype(int)-1])
            tick_text.append(c)
            graph_box_before_normalization.append(plotly_graph_data_before_normalization[c])
            graph_box_after_normalization.append(plotly_graph_data_after_normalization[c])
        f = save_df(pd_from_r_df, "normalization;table")
        graphData = File(columns=[], file_type="plotly;profile_plot;boxplot;json")
        graphData.file.save(f"{str(graphData.link_id)}.json", ContentFile(dumps({
            'tick_val': tick_val, 'tick_text': tick_text, 'graph_box_before_normalization': graph_box_before_normalization,
            'graph_box_after_normalization': graph_box_after_normalization}).encode("utf-8")))
        graphData.save()

        o.output_files.set([f,graphData])
        o.job_finished = True
        o.save()
        message_template["message"] = "Completed operation"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
    except Exception as error:
        o.job_error = str(error)
        o.job_error_status = True
        o.save()
        message_template["message"] = "Error: " + str(error)
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })


def load_dataframe(channel_layer, input_file_id, message_template, request_form, session_id):
    if type(input_file_id) == int:
        file = File.objects.get(id=input_file_id)
    else:
        file = input_file_id
    extension = file.file.name.split(".")[-1]
    async_to_sync(channel_layer.group_send)(session_id, {
        'type': 'job_message',
        'message': message_template
    })
    message_template["message"] = "Loading file"
    async_to_sync(channel_layer.group_send)(session_id, {
        'type': 'job_message',
        'message': message_template
    })
    if extension == "csv":
        df = pd.read_csv(file.file.url)
    elif extension == "xlsx":
        df = pd.read_excel(file.file.url)
    elif extension == "txt" or extension == "tsv":
        df = pd.read_csv(file.file.url, sep="\t")
    df.replace("#VALUE!", np.NAN, inplace=True)
    df[request_form['sampleColumns']] = df[request_form['sampleColumns']].astype(float)
    message_template["message"] = "Convert dataframe to R dataframe"
    async_to_sync(channel_layer.group_send)(session_id, {
        'type': 'job_message',
        'message': message_template
    })
    return df


@job('default')
def r_correlation_matrix(operation_id: int, session_id: str):
    o = Operation.objects.get(id=operation_id)
    message_template = {
        'message': "Running operation",
        'senderName': "Server",
        'requestType': "Correlation Matrix (R)",
        'operationId': o.id
    }
    channel_layer = get_channel_layer()
    request_form = loads(o.value)
    input_file_id = o.input_files.all()[0].id
    try:
        df = load_dataframe(channel_layer, input_file_id, message_template, request_form, session_id)

        with (ro.default_converter + pandas2ri.converter).context():
            df.rename(columns={request_form["indexColumn"]: replace_special_with_dot(request_form["indexColumn"])}, inplace=True)
            r_from_pd_df = ro.conversion.get_conversion().py2rpy(df)
            corrplot = importr("corrplot")

        message_template["message"] = "Creating correlation matrix"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        ro.globalenv["data"] = r_from_pd_df
        ro.globalenv["selectedColumns"] = ro.StrVector(request_form["sampleColumns"] + [replace_special_with_dot(request_form["indexColumn"])])
        ro.r(f"""
        
        data <- data[selectedColumns]
        rownames(data) <- data${replace_special_with_dot(request_form["indexColumn"])}
        data${replace_special_with_dot(request_form["indexColumn"])} <- NULL
        """)

        if request_form["correlationMethod"] in ["pearson", "kendall", "spearman"]:
            ro.r(f"""
            data <- cor(data, method="{request_form["correlationMethod"]}")
            """)

        ro.r("""
        data[is.na(data)] <- 1
        """)

        if request_form["minValue"]:
            ro.r(f"""
            min.value <- {request_form["minValue"]}
            """)
        else:
            ro.r("""
            min.value <- round(min(data), 1) - 0.1
            """)

        ro.r("""
        max.value <- 1
        """)
        if not request_form["correlationMethod"]:
            if request_form["order"] == "hclust":
                ro.r(f"""
                pdf("{input_file_id}_col_corr.pdf")
                cor.data <- corrplot(as.matrix(data), order="{request_form['order']}", hclust.method="{request_form['hclusteringMethod']}", method="{request_form['presentingMethod']}", type="{request_form['correlationPlotShape']}", is.corr=FALSE, col.lim = c(min.value,max.value))
                dev.off()
                """)
            else:
                ro.r(f"""
                pdf("{input_file_id}_col_corr.pdf")
                cor.data <-corrplot(as.matrix(data), order="{request_form['order']}", method="{request_form['presentingMethod']}", type="{request_form['correlationPlotShape']}", is.corr=FALSE, col.lim = c(min.value,max.value))
                dev.off()
                """)
            result = ro.r("cor.data$corrPos")
        else:
            if request_form["order"] == "hclust":
                ro.r(f"""
                pdf("{input_file_id}_col_corr.pdf")
                cor.data <-corrplot(data, order="{request_form['order']}", hclust.method="{request_form['hclusteringMethod']}", method="{request_form['presentingMethod']}", type="{request_form['correlationPlotShape']}", is.corr=FALSE, col.lim = c(min.value,max.value))
                dev.off()
                """)
            else:
                ro.r(f"""
                pdf("{input_file_id}_col_corr.pdf")
                cor.data <- corrplot(data, order="{request_form['order']}", method="{request_form['presentingMethod']}", type="{request_form['correlationPlotShape']}", is.corr=FALSE, col.lim = c(min.value,max.value))
                dev.off()
                """)

            result = ro.r("cor.data$corrPos")
        with (ro.default_converter + pandas2ri.converter).context():
            pd_from_r_df = ro.conversion.get_conversion().rpy2py(result)

        corr_f = save_df(pd_from_r_df, "correlation_matrix;table")
        f = File(columns=[], file_type="correlation_matrix;plot;pdf")
        with open(f"{o.input_files.all()[0].id}_col_corr.pdf", "rb") as pdf_file:
            f.file.save(f"{str(f.link_id)}.pdf", ContentFile(pdf_file.read()))
        f.save()
        o.output_files.set([corr_f, f])
        o.job_finished = True
        o.save()
        message_template["message"] = "Completed operation"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        if os.path.exists(f"{input_file_id}_col_corr.pdf"):
            os.remove(f"{input_file_id}_col_corr.pdf")

    except Exception as error:
        if os.path.exists(f"{o.input_files.all()[0].id}_col_corr.pdf"):
            os.remove(f"{o.input_files.all()[0].id}_col_corr.pdf")

        o.job_error = str(error)
        o.job_error_status = True
        o.save()
        message_template["message"] = "Error: " + str(error)
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })

@job('default')
def convert_msfragger_to_curtainptm(operation_id: int, session_id: str):
    o = Operation.objects.get(id=operation_id)
    message_template = {
        'message': "Running operation",
        'senderName': "Server",
        'requestType': "Correlation Matrix (R)",
        'operationId': o.id
    }
    channel_layer = get_channel_layer()
    request_form = loads(o.value)
    input_file_id = 0
    fasta_file_id = 0
    for i in o.input_files.all():
        if i.file_type == "msfragger;table":
            input_file_id = i.id
        elif i.file_type == "fasta;txt":
            fasta_file_id = i.id

    try:
        df = load_dataframe(channel_layer, input_file_id, message_template, request_form, session_id)
        message_template["message"] = "Parse fasta into dictionary"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        fasta_file = File.objects.get(id=fasta_file_id).file.read().decode("utf-8")
        fasta_dict = read_fasta(fasta_file)
        message_template["message"] = "Processing dataframe"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
        for i, row in df.iterrows():
            match = reg_positon_residue.search(row[request_form["indexColumn"]])
            if match:
                position = int(match.group(2))
                residue = match.group(1)
                position_in_peptide = position
                if row["ProteinID"] in fasta_dict:
                    peptide_seq = row["Peptide"].split(";")[0].upper()
                    try:
                        peptide_position = fasta_dict[row["ProteinID"]].index(peptide_seq)
                    except ValueError:
                        peptide_position = fasta_dict[row["ProteinID"]].replace("I", "L").index(
                            peptide_seq.replace("I", "L"))
                        df.at[i, "Comment"] = "I replaced by L"
                    if peptide_position >= -1:
                        position_in_peptide = position - peptide_position
                df.at[i, "Position"] = position
                df.at[i, "Residue"] = residue
                df.at[i, "Position.in.peptide"] = position_in_peptide

        file_type = "msfragger;curtain;table"
        f = save_df(df, file_type)
        o.output_files.set([f])
        o.job_finished = True
        o.save()
        message_template["message"] = "Completed operation"
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })
    except Exception as error:
        o.job_error = str(error)
        o.job_error_status = True
        o.save()
        message_template["message"] = "Error: " + str(error)
        async_to_sync(channel_layer.group_send)(session_id, {
            'type': 'job_message',
            'message': message_template
        })


def save_df(df, file_type):
    str_data = df.to_csv(sep="\t", index=False)
    f = File(columns=dumps(df.columns.tolist()), file_type=file_type)
    f.file.save(f"{str(f.link_id)}.txt", ContentFile(str_data.encode("utf-8")))
    f.save()
    return f