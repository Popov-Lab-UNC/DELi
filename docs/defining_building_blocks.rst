========================
Defining Building Blocks
========================

Building Blocks are simply just a Building block ID and the corresponding SMILES
::
    {
        "id": "BB_967",
        "smiles": "CCCNCCO"
    }

A Building Block set is just a list of all the building blocks
::
    [
        {
            "id": "BB_967",
            "smiles": "CCCNCCO"
        },
        {
            "id": "BB_145",
            "smiles": "C1CCCCC1"
        },
        {
            "id": "BB_498",
            "smiles": "CCC(=O)O"
        },
    ]

File Format
===========

Building Blocks files should be saved in the the ``$BB_SET_FOLDER`` and have the same name as their id

For example
::
    {
        "id": "BB_A",
        "building_blocks": [
            {
                "id": "BB_967",
                "smiles": "CCCNCCO"
            },
            {
                "id": "BB_145",
                "smiles": "C1CCCCC1"
            },
            {
                "id": "BB_498",
                "smiles": "CCC(=O)O"
            },
        ]
    }

and the file should be named ``BB_A.json``
