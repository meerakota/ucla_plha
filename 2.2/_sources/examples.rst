Examples
========

A set of example configuration files is provided to illustrate the key features of
the **ucla_plha** code, ranging from a single GMM based PSHA to a full PSHA + PLHA calculation
with disaggregation across multiple ground motion and triggering models. 

Each example config is summarized in the table below; the file contents can be
expanded beneath it.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - config file
     - description
   * - `config1.json <https://github.com/meerakota/ucla_plha/blob/main/examples/config1.json>`_
     - PSHA (no PLHA) using BSSA14
   * - `config2.json <https://github.com/meerakota/ucla_plha/blob/main/examples/config2.json>`_
     - PSHA (no PLHA) using ASK14
   * - `config3.json <https://github.com/meerakota/ucla_plha/blob/main/examples/config3.json>`_
     - PSHA (no PLHA) using CB14
   * - `config4.json <https://github.com/meerakota/ucla_plha/blob/main/examples/config4.json>`_
     - PSHA (no PLHA) using CY14
   * - `config5.json <https://github.com/meerakota/ucla_plha/blob/main/examples/config5.json>`_
     - PSHA (no PLHA) using BSSA14, ASK14, CB14, and CY14
   * - `config6.json <https://github.com/meerakota/ucla_plha/blob/main/examples/config6.json>`_
     - PSHA with disaggregation using BSSA14, ASK14, CB14, and CY14
   * - `config7.json <https://github.com/meerakota/ucla_plha/blob/main/examples/config7.json>`_
     - PSHA + PLHA using all four GMMs and the Boulanger & Idriss (2016) triggering model
   * - `config8.json <https://github.com/meerakota/ucla_plha/blob/main/examples/config8.json>`_
     - PSHA + PLHA using all four GMMs and the Moss et al. (2006) triggering model
   * - `config9.json <https://github.com/meerakota/ucla_plha/blob/main/examples/config9.json>`_
     - PSHA + PLHA using all four GMMs and the NGL SMT triggering model
   * - `config10.json <https://github.com/meerakota/ucla_plha/blob/main/examples/config10.json>`_
     - PSHA + PLHA with disaggregation, all four GMMs, and three triggering models (Boulanger & Idriss 2016, NGL SMT, Moss et al. 2006); also writes an output file

Configuration file contents
---------------------------

.. dropdown:: config1.json -- PSHA (no PLHA) using BSSA14

   .. literalinclude:: ../examples/config1.json
      :language: json

.. dropdown:: config2.json -- PSHA (no PLHA) using ASK14

   .. literalinclude:: ../examples/config2.json
      :language: json

.. dropdown:: config3.json -- PSHA (no PLHA) using CB14

   .. literalinclude:: ../examples/config3.json
      :language: json

.. dropdown:: config4.json -- PSHA (no PLHA) using CY14

   .. literalinclude:: ../examples/config4.json
      :language: json

.. dropdown:: config5.json -- PSHA (no PLHA) using BSSA14, ASK14, CB14, and CY14

   .. literalinclude:: ../examples/config5.json
      :language: json

.. dropdown:: config6.json -- PSHA with disaggregation using BSSA14, ASK14, CB14, and CY14

   .. literalinclude:: ../examples/config6.json
      :language: json

.. dropdown:: config7.json -- PSHA + PLHA using all four GMMs and the Boulanger & Idriss (2016) triggering model

   .. literalinclude:: ../examples/config7.json
      :language: json

.. dropdown:: config8.json -- PSHA + PLHA using all four GMMs and the Moss et al. (2006) triggering model

   .. literalinclude:: ../examples/config8.json
      :language: json

.. dropdown:: config9.json -- PSHA + PLHA using all four GMMs and the NGL SMT triggering model

   .. literalinclude:: ../examples/config9.json
      :language: json

.. dropdown:: config10.json -- PSHA + PLHA with disaggregation, all four GMMs, and three triggering models (Boulanger & Idriss 2016, NGL SMT, Moss et al. 2006); also writes an output file

   .. literalinclude:: ../examples/config10.json
      :language: json

Worked notebook
---------------

A complete, runnable walkthrough is available as a Jupyter notebook,
`Example_ucla_plha_notebook.ipynb <https://github.com/meerakota/ucla_plha/blob/main/examples/Example_ucla_plha_notebook.ipynb>`_, which
loads a configuration, runs the analysis, and plots hazard and disaggregation results.
