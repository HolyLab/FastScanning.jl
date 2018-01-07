function fake_imgs(nconditions, ntrials, best_condition; best_img::AbstractArray{T,2} = rand(2,2)) where {T}
    output = rand(size(best_img)..., nconditions*ntrials)
    fake_imgs!(output, nconditions, ntrials, best_condition; best_img = best_img)
end

function fake_imgs!(output_img::AbstractArray{T,3}, nconditions, ntrials, best_condition; best_img::AbstractArray{T,2} = rand(2,2)) where {T}
    @assert size(output_img)[1:2] == size(best_img)
    nimgs = nconditions * ntrials
    @assert isinteger(nconditions)
    for c = 1:nconditions
        if c == best_condition
            output_img[:,:, Int((c-1)*ntrials+1):2:Int((c*ntrials))] = best_img #at optimal lag
        else
            output_img[:,:, Int((c-1)*ntrials+1):2:Int((c*ntrials))] = rand(size(best_img)...)
        end
    end
    return output_img
end

function fake_slicetiming_run(nconditions, ntrials, nslices, best_conds_fwd, best_conds_back)
    @assert nslices == length(best_conds_fwd) == length(best_conds_back)
    output = zeros(2,2,2*nconditions*ntrials*nslices + nslices)
    output_fwd = view(output, :,:,nslices+1:2:size(output,3))
    output_back = view(output, :,:,nslices+2:2:size(output,3))
    template =rand(2,2,nslices)
    output[:,:,1:nslices] = template
    for i = 1:nslices
        fwd_sub = view(output_fwd, :,:,1+(i-1)*nconditions*ntrials:i*nconditions*ntrials)
        fake_imgs!(fwd_sub, nconditions, ntrials, best_conds_fwd[i]; best_img = template[:,:,i])
        back_sub = view(output_back, :,:,1+(i-1)*nconditions*ntrials:i*nconditions*ntrials)
        fake_imgs!(back_sub, nconditions, ntrials, best_conds_back[i]; best_img = template[:,:,i])
    end
    return output
end
