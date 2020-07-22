Reproducing the Acceptance Report
=================================

1. Download the document as a .docx file: https://docs.google.com/document/d/158ba_xQM9AYyu8VcLWsIvSoEYps6PQhgddTr9H0BFmY/edit?ts=5c876d41#
2. Ensure you have pandoc installed on your system.
3. Run the following code to convert the results to .rst (this code must be run on a POSIX system).

.. code-block:: bash

    # convert the file
    pandoc -f docx Optimus\ Acceptance\ Report.docx -t rst -o optimus_report.rst

    # get the media (image) files from the .docx file
    unzip ./Optimus\ Acceptance\ Report.docx
    mv ./word/media media
