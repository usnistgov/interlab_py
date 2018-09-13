Subobjects within the Project class
***********************************

This page presents documentation for the subobjects that are contained within a :py:class:`.Project` object. All of the methods of these objects can be accessed from the :py:class:`.Project` object, so the documentation here is sparse. The :py:class:`.Project` object will contain :py:class:`.ExperimentGroup` and :py:class:`.InterlabArray` objects. :py:class:`.ExperimentGroup` objects will, in turn, contain :py:class:`.DistanceMetric` objects, and both :py:class:`.InterlabArray` and :py:class:`.DistanceMetric` objects will contain :py:class:`.Population` objects. Container objects usually contain an interface to the equivalent method in the contained object.

.. currentmodule:: utilities

Class documentation
-------------------
   
.. autoclass:: ExperimentGroup
.. autoclass:: InterlabArray
.. autoclass:: DistanceMetric
.. autoclass:: Population

Method summary
--------------

.. autosummary::
   ExperimentGroup
   ExperimentGroup.distance
   ExperimentGroup.fit_zscores
   ExperimentGroup.find_outliers
   ExperimentGroup.histogram
   ExperimentGroup.plot_data
   ExperimentGroup.plot_zscores
   ExperimentGroup.distance_measure_plot
   
.. autosummary::
   InterlabArray
   InterlabArray.fit_transform
   InterlabArray.fit_zscores
   InterlabArray.find_outliers
   InterlabArray.plot_components
   InterlabArray.plot_zscores
   InterlabArray.plot_outliers

.. autosummary::
   DistanceMetric
   DistanceMetric.distance
   DistanceMetric.fit_zscores
   DistanceMetric.find_outliers
   DistanceMetric.plot_zscores
   
.. autosummary::
   Population
   Population.fit_zscores
   Population.find_outliers
   Population.histogram
   Population.plot_zscores
   
   