FROM python:3.10-bullseye
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
LABEL authors="Toan Phung"

RUN curl https://www.postgresql.org/media/keys/ACCC4CF8.asc | gpg --dearmor | tee /etc/apt/trusted.gpg.d/apt.postgresql.org.gpg > /dev/null
RUN sh -c 'echo "deb http://apt.postgresql.org/pub/repos/apt bullseye-pgdg main" > /etc/apt/sources.list.d/pgdg.list'
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 7FCC7D46ACCC4CF8
RUN apt-get update
RUN apt-get -y install postgresql-client-14
WORKDIR /app
RUN mkdir "/app/media"
RUN mkdir "/app/backup"
RUN apt-get install -y --no-install-recommends r-base r-base-dev && rm -rf /var/lib/apt/lists/*
RUN R -e 'install.packages("devtools")'
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('remotes', 'dplyr', 'ggplot2'))"
RUN R -e "library('BiocManager');BiocManager::install(c('MSstats', 'MSstatsPTM', 'limma', 'preprocessCore', 'impute', 'QFeatures', 'norm', 'pcaMethods', 'imputeLCMD', 'corrplot', 'PerseusR', 'data.table', 'XML', 'doParallel', 'stringr', 'MASS', 'Rcpp', 'Biostrings', 'Biobase'))"
COPY . /app/
RUN pip install -r requirements.txt
EXPOSE 8000

CMD ["gunicorn", "currantDjango.wsgi:application", "--timeout", "300", "--bind", "0.0.0.0:8000"]