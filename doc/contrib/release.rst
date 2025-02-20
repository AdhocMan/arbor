Releases
********

Release cycle
=============

0. release every 3 months (at time ``T``)
1. ``T-11`` weeks: ``all`` add your favorite Issues to the next-rel column
2. ``T-10`` weeks: ``Scrum Master`` prep dev meet (internal)

   * Update/trim next-release column in Kanban
   * Prepare agenda, include possible additions not covered by Kanban/Issues
   * Add milestone tags (nextver, nextver+1, etc.)
   * Add highlighted new features to RELEASE_NOTES.md
3. ``T-8`` weeks: ``Release Manager`` dev meet (external/public)

   * Use Kanban as starter
   * Move issues around based on input
   * Add milestone tags, for this release or future releases
   * Update highlighted new features to RELEASE_NOTES.md
4. ``T-1``: ``all`` reserve week for wrapping up PRs and review.
5. ``T±0``: ``Release Manager`` release!

   * Have a look at Python release schedule, and time Arbor release optimally with new Python minor version. It is nice to generate wheels for the new minor as soon as minor is released.
6. ``T+1`` weeks: ``Scrum Master`` retrospective
   
   * set date for next release

Procedure
=========

These notes enumerate the steps required every time we release a new
version of Arbor.

Pre-release
-----------

Update tags/versions and test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Check if some files are up to date
    
   - README.md, ATTRIBUTIONS.md, CONTRIBUTING.md, RELEASE_NOTES.md
   - Verify MANIFEST.in (required for PyPI sdist)
   - Double check that all examples/tutorials/etc are covered by CI
   - Check Python/pip/PyPi metadata and scripts, e.g. ``setup.py``, ``pyproject.toml``

#. Create new temp-branch ending in ``-rc``. E.g. ``v0.6-rc``
#. Bump the ``VERSION`` file:

   - See also :ref:`dev-version`
   - Append ``-rc``. (Make sure there's no ``-dev``)

#. Run all tests.

   - ``ciwheel.yml`` triggers when you push a branch called ``ciwheel``, and on new Git tags. Since we've not tagged the release yet, run ``git push origin HEAD:ciwheel``. Make sure the tests pass.
   
      - ``ciwheel.yml`` pushes automatically to `Test.PyPI.org <https://test.pypi.org/project/arbor/>`_. This only passes if ran as branch of the main ``arbor-sim`` repo (as that's where the PyPI secret lives). On your own repo, the upload will fail (the rest should pass). If you want to test uploading, then force push to the _upstream_ ``ciwheel`` branch, e.g. ``git push upstream HEAD:ciwheel --force``.
   
   - If you want to test the PyPI submission (consider asking other OS-users):

     .. code-block:: bash

        python -m venv env && source env/bin/activate
        pip install numpy
        pip install -i https://test.pypi.org/simple/ arbor #should select the latest build, but doublecheck
        python -c 'import arbor; print(arbor.__config__)'

   - Use build flags to test the source package: :ref:`in_python_adv`

Release
-------

#. Make sure no errors were encountered in the pre-release phase, working wheels were produced, and were installable from `Test.PyPI.org <https://test.pypi.org/project/arbor/>`_ and that no problems were reported.
   
#. Change ``VERSION``. Make sure does not end with ``-rc`` or ``-dev``.

#. Update ``scripts/check-all-tags.sh`` to check the current tag.

#. Tag

   - commit and open a PR for the above changes.
   - on cmdline: ``git tag -a TAGNAME``
   - ``git push upstream TAGNAME``

#. Upload to pypi & verify

   .. code-block:: bash

      twine upload -r arborpypi dist/*

      python -m venv env && source env/bin/activate
      pip install arbor
      python -c 'import arbor; print(arbor.__config__)'

#. Create tarball with
   ``scripts/create_tarball ~/loc/of/arbor tagname outputfile``

   - eg ``scripts/create_tarball /full/path/to/arbor v0.5.1 ~/arbor-v0.5.1-full.tar.gz``
   
#. Download output of wheel action associated to this release commit and extract (verify the wheels and
   source targz is in /dist)

   - Of course, the above action must have passed the tests successfully.
   
#. Update ``spack/package.py``. The checksum of the targz is the sha256sum.

#. Start a new release on Zenodo, this allocated a DOI, but you don't have to finish it right away. Add new Zenodo badge/link to docs/README.

#. Create Github Release: https://github.com/arbor-sim/arbor/releases

   - Go to `GH tags`_ and click “…” and “Create release”
   - Categorize/edit Github's autogenerated release notes (alternatively go through merged PRs to come up with a changelog).
   - add tarball to release, created in previous step.
   
#. Update Zenodo with authors and changelog created in previous step and submit.

Post Release
------------

#. Make a new PR setting ``VERSION`` to the next with a trailing ``-dev``. E.g. if you just release ``3.14``, change ``VERSION`` to ``3.15-dev``
    
   - Include changes such as to ``spack/package.py``, ``CITATIONS``, ``doc/index.rst`` in postrel PR. Copy Zenodo BibTex export to ``CITATIONS``.

#. Update spack package / Ebrains Lab / Opensourcebrain

   - Spack upstream: `PR here <https://github.com/spack/spack/blob/develop/var/spack/repos/builtin/packages/arbor/package.py>`_
   - Ebrains Lab: `MR here <https://gitlab.ebrains.eu/technical-coordination/project-internal/devops/platform/ebrains-spack-builds/>`_
   - OSB: update `dockerfile <https://github.com/OpenSourceBrain/OSBv2/blob/master/applications/jupyterlab/Dockerfile>`_ if needed.

     - Make sure that `Notebooks <https://www.v2.opensourcebrain.org/repositories/38>`_ work on the version that their image is built with.

#. Announce on our website
#. Announce on HBP newsletter newsletter@humanbrainproject.eu, HBP Twitter/socials evan.hancock@ebrains.eu
#. [AUTOMATED] Add tagged version of docs on ReadTheDocs
#. HBP internal admin

   - Plus: https://plus.humanbrainproject.eu/components/2691/
   - TC Wiki: https://wiki.ebrains.eu/bin/view/Collabs/technical-coordination/EBRAINS%20components/Arbor/
   - KG: https://search.kg.ebrains.eu/instances/5cf4e24b-b0eb-4d05-96e5-a7751134a061
 
     - Update howto: https://github.com/bweyers/HBPVisCatalogue/wiki/How-to-start-software-meta-data-curation%3F#update-curated-software
     - Previous update as template: https://github.com/bweyers/HBPVisCatalogue/issues/480
     - Supported file formats
 
       - ContentTypes: https://humanbrainproject.github.io/openMINDS/v3/core/v4/data/contentType.html
       - details: https://github.com/HumanBrainProject/openMINDS_core/tree/v3/instances/data/contentTypes
 
   - Send an update to the folk in charge of HBP Twitter if we want to shout about it

#. FZJ admin

   - https://juser.fz-juelich.de/submit

.. _GH tags: https://github.com/arbor-sim/arbor/tags
.. _AUTOMATED: https://github.com/arbor-sim/arbor/blob/master/.github/workflows/ebrains.yml 
