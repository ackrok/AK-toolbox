function writeNeuroscopeShank(clusters,samplingRate,fileBaseName)
% takes cluster IDs and spike times from loadKS (or variant) and writes two file types:
% fileBaseName.clu.n = list of cluster IDs in chronological order
%                       corresponding to spike times; first entry = total
%                       number of clusters on shank (n = shank number)
% fileBaseName.res.n = list of spike times in chronological order;
%                       (n = shank number)

numClust = length(clusters);
shankID = horzcat(clusters(:).shank);
clusterID = double(horzcat(clusters(:).clusterID));
numShank = unique(shankID);


for n = numShank
    indN = find(shankID == n);
    shankSpikeTimes = [];
    shankClusterID = [];
    spikeTimes = [];
    clustIDvec = [];
    shankClusterIDSorted = [];
    shankSpikeTimesSorted = [];
    for k = indN
        spikeTimes = clusters(k).spikeTimes.*samplingRate;
        clustIDvec = ones(length(spikeTimes),1);
        clustIDvec = clustIDvec.*clusterID(k);
        shankSpikeTimes = [shankSpikeTimes;spikeTimes]; 
        shankClusterID = [shankClusterID;clustIDvec];
    end
    [shankSpikeTimesSorted,indsT] = sort(shankSpikeTimes);
    shankClusterIDSorted = shankClusterID(indsT);    
    shankClusterIDSorted = [length(indN);shankClusterIDSorted];
    dlmwrite([fileBaseName,'.clu.',num2str(n)],shankClusterIDSorted,'precision',11,'newline','pc');
    dlmwrite([fileBaseName,'.res.',num2str(n)],shankSpikeTimesSorted,'precision',11,'newline','pc');
    
end

end

