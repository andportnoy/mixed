DMATRICES = X Z Lambdat
FMATRICES = Whalf WX Wy ZtW XtWX XtWy ZtWX ZtWy L

DMATRICES_R  = $(addsuffix -r.bin,  $(DMATRICES))
DMATRICES_PY = $(addsuffix -py.bin, $(DMATRICES))

FMATRICES_R  = $(addsuffix -r.bin,  $(FMATRICES))
FMATRICES_PY = $(addsuffix -py.bin, $(FMATRICES))

# we can get -0 vs +0, so binary comparison won't work here
define cmp
@echo 'Comparing $1 in Python and R...'
@python3 cmp.py -a $1-py.bin -b $1-r.bin

endef

check: check_design check_fit

check_fit: lme4pureR $(FMATRICES_R) $(FMATRICES_PY)
	$(foreach mat,$(FMATRICES),$(call cmp,$(mat)))

$(FMATRICES_PY): data.feather formula.txt
	@echo 'Generating fitting matrices in Python...'
	@python3 run_pls.py --formula formula.txt --data data.feather

$(FMATRICES_R): data.feather formula.txt
	@echo 'Generating fitting matrices in R...'
	@Rscript run_pls.r --formula formula.txt --data data.feather 2> /dev/null

.PHONY: lme4pureR
lme4pureR:
	@echo 'Reinstalling lme4pureR...'
	@R CMD INSTALL --no-help --no-byte-compile --no-test-load $@ 2> /dev/null

check_design: $(DMATRICES_PY) $(DMATRICES_R)
	$(foreach mat,$(DMATRICES),$(call cmp,$(mat)))

$(DMATRICES_R): data.feather formula.txt
	@echo 'Generating design matrices in R...'
	@Rscript build_matrices.r --data data.feather --formula formula.txt \
		$(foreach m,$(DMATRICES),--$m $m-r.bin)

$(DMATRICES_PY): data.feather formula.txt
	@echo 'Generating design matrices in Python...'
	@python3 build_matrices.py --data data.feather --formula formula.txt \
		$(foreach m,$(DMATRICES),--$m $m-py.bin)

data.feather: columns.txt
	@echo 'Generating data from list of columns...'
	@python3 generate_data.py --columns $^ --data $@

setup:
	pip3 install -U patsy pandas feather-format
	MAKE='make -j' Rscript -e 'install.packages(c("lme4", "optparse", "feather"), repos="cloud.r-project.org")'

clean:
	rm -rf $(wildcard *.feather *.bin __pycache__)
