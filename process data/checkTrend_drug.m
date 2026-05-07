% check on signal
%%
win = 15*60;
nFPchan = 2;
opts = {'exp2stitch','exp2base','stitch','stretch'};

for b = 1:nFPchan
    for a = 1:length(comb)
        if isfield(comb,'nbFP'); signal = comb(a).nbFP{b};
            else, signal = comb(a).FP{b}; 
        end
        Fs = comb(a).Fs;
        out = detrend_drug(signal, Fs, win);
        [ax, fig] = plotTrend_drug(out);
        ID = [comb(a).mouse,'-',comb(a).date];
        choice = menu(sprintf('Select trend for dF/F: %s',ID),opts);
        chosen = opts{choice};
        dff = getdff_drug(out.y, out.trend.(chosen), Fs, win);
        comb(a).dff(:,b) = dff;
        comb(a).trend{b} = chosen; 
        close(fig);
    end
end