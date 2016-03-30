function out = logit(x)
    out=(-log(1./x-1));
    out(~isfinite(out))=NaN;
end