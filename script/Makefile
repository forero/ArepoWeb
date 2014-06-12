IN_PATH=/home/springel/universe/AnalysisIllustris/AddGrads/L75n455FP/snapdir_136/
OUT_PATH=/home/extforer/ArepoWeb/data/
SNAPS = $(notdir $(wildcard $(IN_PATH)snap_136*.hdf5))
TARGETS = $(subst snap,gasweb_snap,$(SNAPS))
EXEC=/home/extforer/ArepoWeb/src/./snap2web.x

all: $(TARGETS)

$(TARGETS): gasweb_%: $(IN_PATH)%
	ls -l $<
	touch $@
	$(EXEC) $< $@

clean:
	rm -f $(TARGETS)