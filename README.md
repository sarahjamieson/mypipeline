# Myeloid-Pipeline-STP

A Django project for running an NGS pipeline and storing the data in a database, then presenting the result in a web app. This project is specific to MiSeq Nextera data.

The pipeline detects common structural rearrangements and FLT3 internal tandem duplications (FLT3-ITDs) associated with acute myeloid leukaemia (AML).

The pipeline takes 3 arguments:

| Argument      | Description                                                                 |
|---------------|-----------------------------------------------------------------------------|
|-s/--sheet     | A path to the Illumina sample sheet for this run (CSV file).                |
|-d/--result_dir| A path to the directory containing the MiSeq InterOp and Data directories.  |
|-o/--output_dir| A path for the results to be stored.                                        |


To run the pipeline:
```
$ python pipeline.py -s /path/to/samplesheet.csv -d /path/to/miseq/output/ -o /path/to/result/
```

To deploy the app:

1. cd to the directory with the `manage.py` file.
2. `$ python manage.py runserver 127.0.0.1:8000`

`127.0.0.1:8000` is a development server; to deploy on another server replace this with another IP address.
  
