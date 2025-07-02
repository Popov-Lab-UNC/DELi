.. _selection-file-docs:

===============
Selection Files
===============
DEL selections are experiments, and like all experiments, they have all kinds of metadata and conditions that are associated with them.
Some of this information is useful for DELi to know when processing the selection data. To help streamline this DELi has a "selection file
format" that allows you to easily enter information about your DEL selection. It is based on YAML and looks like this:

.. code-block:: yaml

    selection_id: "EXAMPLE_SELECTION"
    target_id: "Target1"
    selection_condition: "200uL Library with 50pM Protein"
    date_ran: "2025-03-30"
    additional_info: "This is an example selection for DELi documentation."

    libraries:
        - "DEL004"
        - "DEL005"

Below is a description of the fields in this file:

- **selection_id**: An identifier for the selection. This should be a short, descriptive name that identifies the selection
  (preferably the same as any other experiment tracking you use).
- **target_id**: The identifier for the target protein or molecule that the selection is being performed against.
  This should match the target ID used in any of your other traget tracking, or could be something like the Uniprot ID.
- **selection_condition**: A description of the conditions under which the selection was performed.
  This could include information about the library concentration, the target protein concentration, or any other relevant details.
- **date_ran**: The date the selection was performed. This should be in ISO 8601 format (YYYY-MM-DD).
- **additional_info**: Any additional information about the selection that you want to include.
  This could include notes about the selection process, any issues that were encountered, or any other relevant details worth logging.

After this metadata/experiment information, you then list the libraries that were used in the selection. DELi will use the library
information to help process the selection data, both in a the analysis module and the decoding module.

.. note::
    You can use the full path to the :ref:`library JSON files <defining_libraries>`, or just the name of the library if it is in
    the :ref:`DELi Data Directory <deli-data-dir-ref>`.

These selection files can also serve as a good way to track your experiments. While the list keys above are the only required ones,
(except for ``additional_info``, which is optional), you can add any other keys you want to track to the file.
DELi will ignore any keys that it does not recognize, so you can add as many as you want.
