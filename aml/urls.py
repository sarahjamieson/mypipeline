from django.conf.urls import url
from django.contrib import admin
from aml.views import index, get_samples_for_run, get_results_for_sample

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^$', index, name='index'),
    url(r'^(?P<run>[-\w]+)/$', get_samples_for_run, name='samples'),
    url(r'^(?P<sample>[-\w]+)/$', get_results_for_sample, name='results'),
]
