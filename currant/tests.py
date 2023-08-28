from django.test import TestCase
from django_rq import get_worker

import os
os.environ["R_HOME"] = "C:\\Program Files\\R\\R-4.3.1"

from currant.tasks import save_file, r_qfeatures_protein
from currantDjango.settings import RQ_QUEUES


class TestRQ(TestCase):
    def test_rq(self):

        for queueConfig in RQ_QUEUES.values():
            queueConfig['ASYNC'] = False
        request_form = {
            "normalization": "quantiles.robust",
            "imputation": "knn",
            "dataCompleteness": "1",
            "indexColumn": "T: Index",
            'sampleColumns': ["4Hr-AGB1.01","4Hr-AGB1.02","4Hr-AGB1.03","4Hr-AGB1.04","4Hr-AGB1.05","24Hr-AGB1.01","24Hr-AGB1.02","24Hr-AGB1.03","24Hr-AGB1.04","24Hr-AGB1.05","4Hr-Cis.01","4Hr-Cis.02","4Hr-Cis.03","24Hr-Cis.01","24Hr-Cis.02","24Hr-Cis.03"],
            "conditionMap": {},
            "log2": True
        }

        for s in request_form['sampleColumns']:
            request_form['conditionMap'][s] = {"condition": s.split(".")[0], "replicate": s.split(".")[1]}
        r_qfeatures_protein.delay(19, request_form)
        worker = get_worker()
        worker.work(burst=True)
        self.assertTrue(True)
