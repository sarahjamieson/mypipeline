import os
import datetime

sens = [91, 90]
log_file = open("multiple_log.txt", "w")  # write into output file

for item in sens:
    start = datetime.datetime.now().isoformat()
    log_file.write("%s: Starting sensitivity 0.%s." % (start, item))
    os.system("python /home/shjn/PycharmProjects/mypipeline/aml/pipeline.py "
              "-s /media/sf_sarah_share/160801_s%s/SampleSheet.csv "
              "-d /media/sf_sarah_share/160801_s%s/ "
              "-o /media/sf_sarah_share/MiSeq_Nextera_Results/ "
              "-p 0.%s"
              % (item, item, item))
    end = datetime.datetime.now().isoformat()
    log_file.write("%s: Ending sensitivity %s." % (end, item))

log_file.close()
