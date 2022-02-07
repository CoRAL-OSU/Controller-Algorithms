function rsme=rms_error(error)
subs=error.^2;
add=sum(subs);
rsme=sqrt(add/numel(error));
end