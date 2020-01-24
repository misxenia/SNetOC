function [samples] = thinning(samples, thin)

names = fieldnames(samples);
nsamples = size(samples(1).(names{1}), ndims(samples(1).(names{1})));
[ntypes, nchains] = size(samples);

ind = thin:thin:nsamples;

for t=ntypes:-1:1 %% loop over types of nodes
    for j=1:numel(names)
        name = names{j};
        if isnumeric(samples(t,1).(name))
            for k=nchains:-1:1
                if ismatrix(samples(t,k).(name))
                    samples(t,k).(name) = samples(t,k).(name)(:, ind);
                else
                    samples(t,k).(name) = samples(t,k).(name)(:, :, ind);
                end
            end
        else
            fn = fieldnames(samples(t,1).(name));
            for v=1:numel(fn)
                for k=nchains:-1:1
                    if ismatrix(samples(t,k).(name).(fn{v}))
                        samples(t,k).(name).(fn{v}) = ...
                            samples(t,k).(name).(fn{v})(:, ind);
                    else
                        samples(t,k).(name).(fn{v}) = ...
                            samples(t,k).(name).(fn{v})(:, :, ind);
                    end
                end
            end
        end
    end
end
    