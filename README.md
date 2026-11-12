# GTHM26


> [!important]
> documentation in progress



The tool can be executed by running the following example calculating hazard curves for the DART 21416 and using the Kurils Japan subduction zone source area:

```
python -m gthm26 --subd kurilsjapan --poi DART-21416 --subd_id 009 --output_path path_to_output_folder --rupt_file path_to_ruptures_file.nc --gf_path path_to_greens_functions_folder --plot
```

where `path_to_output_folder`, `path_to_ruptures_file.nc` and `path_to_greens_functions_folder` depends on where the user place these files in her/his own machine.

The set of sources used in the GTHM26 can be downloaded from:

https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/catalog.html

The user can select a subduction `subname` then click on `TSUNAMI_EVENTS` folder and download the  `all_stochastic_slip_earthquake_events_subname.nc` files in that folder, where `subname`, is the name of the subduction.


Green Functions database can be downloaded from here etc etc etc




