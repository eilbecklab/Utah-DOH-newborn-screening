## LOVD2 Variant Parsing

This is currently a jupyter notebook, but will soon be updated to a python script with command line options that will be read in using the [argparse package](https://docs.python.org/3/library/argparse.html).

Currently, the Cincinnati Children's Hospital Medical Center (CCHMC) LOVD2 database is the only one that we are parsing. The LOVD2 format required visiting the page of each variant to obtain pathogenicity information, so this program goes through the pages for each variant to obtain the information. This will take in lists of genes and scrape the information about the variants within those genes.

The variants have transcript information listed, but the transcripts are often older versions of the transcripts. Biocommons doesn't recognize the older versions of the transcripts, so the provided csv file "Transcript_Info_For_Dictionaries.csv" is provided to allow the program to look up the current version of the transcript for GRCh37. The transcripts tested had no difference between the older and newer transcript versions aside from extending the untranslated regions (UTRs) to further away from the coding regions, indicating that anything between the transcription start site (TSS) and transcription termination site (TTS) will be identical and only the variants far into the UTRs will not be found.
