bin_SCRIPTS = lin_from_nlin scale_voxels add_noise_to_volume grid_object_manipulator ray_trace_crosshair concatenate_tags make_xfm_for_grid.pl

grid_object_manipulator: grid_object_manipulator
	sed < $< > $@-t \
		-e 's,xPKGDATADIRx,$(pkgdatadir),' &&\
	mv $@-t $@

ray_trace_crosshair: ray_trace_crosshair
	sed < $< > $@-t \
		-e 's,xPKGDATADIRx,$(pkgdatadir),' &&\
	mv $@-t $@

EXTRA_DIST = $(bin_SCRIPTS) \
	grid_object_manipulator \
	ray_trace_crosshair
