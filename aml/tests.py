from django.test import TestCase, Client
import django
import os
from django.core.urlresolvers import reverse
import unittest
from selenium import webdriver
os.environ['DJANGO_SETTINGS_MODULE'] = 'mypipeline.settings'
django.setup()


class TestViews(TestCase):
    def setUp(self):
        self.client = Client()

    def test_index(self):
        response = self.client.get(reverse('index'))
        self.assertEqual(response.status_code, 200)
        self.assertIn('<title>Runs</title>', response.content)

    def test_bamstats(self):
        response = self.client.get('/aml/16053/D03-21521/bamstats/')
        self.assertEqual(response.status_code, 200)
        self.assertIn('<title>BAMStats Report</title>', response.content)

    def test_combined(self):
        response = self.client.get(reverse('combined'))
        self.assertEqual(response.status_code, 200)
        self.assertIn('<title>Combined FLT3 Results</title>', response.content)

    def test_delly(self):
        response = self.client.get('/aml/delly/16053/D03-21521/')
        self.assertEqual(response.status_code, 200)
        self.assertIn('<title>Results</title>', response.content)

    def test_fastqc(self):
        response = self.client.get('/aml/16053/D03-21521/fastqc/')
        self.assertEqual(response.status_code, 200)
        self.assertIn('<title>FASTQC Report</title>', response.content)

    def test_flt3(self):
        response = self.client.get('/aml/16053/D03-21521/FLT3/')
        self.assertEqual(response.status_code, 200)
        self.assertIn('<title>FLT3 only results</title>', response.content)

    def test_results(self):
        response = self.client.get('/aml/results/16053/D03-21521/')
        self.assertEqual(response.status_code, 200)
        self.assertIn('<title>Results</title>', response.content)

    def test_samples(self):
        response = self.client.get('/aml/16053/')
        self.assertEqual(response.status_code, 200)
        self.assertIn('<title>16053</title>', response.content)


class TestIntegration(unittest.TestCase):
    """Ensure app deployed on test server 127.0.0.1:8000."""
    def setUp(self):
        self.driver = webdriver.Firefox()

    def test_index(self):
        self.driver.get("http://127.0.0.1:8000/aml/")
        self.assertIn("http://127.0.0.1:8000/", self.driver.current_url)
        self.driver.find_element_by_link_text("16053").click()  # works
        self.assertTrue(self.driver.find_element_by_tag_name("title"), "16053")

    def test_results(self):
        self.driver.get("http://127.0.0.1:8000/aml/results/16053/D15-18331/")
        self.assertIn("http://127.0.0.1:8000/", self.driver.current_url)
        total_rows = self.driver.find_elements_by_xpath("//tr")
        self.assertEqual(len(total_rows), 18)
        element = self.driver.find_element_by_name("gene_filter")
        element.send_keys("FLT3")
        self.driver.find_element_by_xpath("//input[@type='submit']").click()
        filtered_rows = self.driver.find_elements_by_xpath("//tr")
        self.assertEqual(len(filtered_rows), 11)

    def tearDown(self):
        self.driver.quit()

if __name__ == '__main__':
    unittest.main()

# https://realpython.com/blog/python/testing-in-django-part-1-best-practices-and-examples/
# http://selenium-python.readthedocs.io/locating-elements.html
