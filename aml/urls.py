from django.conf.urls import url
from django.contrib import admin
from aml.views import index, get_samples_for_run, get_results_for_sample, get_interop_for_run, get_delly_for_sample, \
    get_sample_quality, get_fastqc, get_bamstats, get_flt3_only

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^$', index, name='index'),
    url(r'^(?P<run>[-\w]+)/$', get_samples_for_run, name='samples'),
    url(r'^results/(?P<run>[-\w]+)/(?P<sample>[-\w]+)/$', get_results_for_sample, name='results'),
    url(r'^delly/(?P<run>[-\w]+)/(?P<sample>[-\w]+)/$', get_delly_for_sample, name='delly'),
    url(r'^interop/(?P<run>[-\w]+)/$', get_interop_for_run, name='interop'),
    url(r'^quality/(?P<run>[-\w]+)/(?P<sample>[-\w]+)/$', get_sample_quality, name='quality'),
    url(r'^(?P<run>[-\w]+)/(?P<sample>[-\w]+)/fastqc/$', get_fastqc, name='fastqc'),
    url(r'^(?P<run>[-\w]+)/(?P<sample>[-\w]+)/bamstats/$', get_bamstats, name='bamstats'),
    url(r'^(?P<run>[-\w]+)/(?P<sample>[-\w]+)/FLT3/$', get_flt3_only, name='flt3'),
]
