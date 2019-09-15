##This is like the toffset_fits function but doesn't evaluate all toffsets:
##does a recursive binary-ish search.  If you want to eliminate a higher or lower fraction of points per-iteration you can adjust elim_frac higher or lower.
##Also differsf from toffset_fits in that it only returns the tfm and mm of the best fit
#function bsearch_toffsets(target, testimgs, toffsets; elim_frac = 0.4, sigmas=(1.0,1.0))
#    @show ntoffsets = length(toffsets)
#    ntrials = convert(Int, size(testimgs,3) / ntoffsets)
#    @show ntrials
#    firstmov = squeeze(mean(Float64.(view(testimgs, :, :, 1:ntrials)), 3),3)
#    firsttfm, firstmm = align2d(target, firstmov; sigmas=sigmas)
#    lastmov = squeeze(mean(Float64.(view(testimgs, :, :, (ntoffsets-1)*ntrials+1:ntoffsets*ntrials)), 3),3)
#    lasttfm, lastmm = align2d(target, lastmov; sigmas=sigmas)
#	if firstmm <= lastmm
#		toffsets_new = toffsets[1: round(Int, ntoffsets*(1-elim_frac))]
#		@show ntoffsets_new =length(toffsets_new)
#		n_elim = ntoffsets - ntoffsets_new
#        @show toffsets_new
#		if ntoffsets_new == 1 || toffsets == toffsets_new
#            return first(toffsets), firsttfm, firstmm
#		else
#    		testimgs = Float64.(view(testimgs, :, :, 1:ntoffsets_new*ntrials))
#            print("here1\n")
#			return bsearch_toffsets(target, testimgs, toffsets_new; elim_frac = elim_frac, sigmas = sigmas)
#		end
#	else
#        toffsets_new = toffsets[(round(Int, ntoffsets*elim_frac)+1):ntoffsets]
#		@show ntoffsets_new = length(toffsets_new)
#		n_elim = ntoffsets - ntoffsets_new
#        @show toffsets_new
#		if ntoffsets_new == 1 || toffsets == toffsets_new
#            return last(toffsets), lasttfm, lastmm
#		else
#            @show n_elim
#            @show ntoffsets
#            @show size(testimgs)
#    		testimgs = Float64.(view(testimgs, :, :, n_elim*ntrials+1:ntoffsets*ntrials))
#            print("here2\n")
#            @show size(testimgs)
#            @show length(toffsets_new)
#			return bsearch_toffsets(target, testimgs, toffsets_new; elim_frac = elim_frac, sigmas = sigmas)
#		end
#	end
#end
#
