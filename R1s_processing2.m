
workDir = '/illumina1/YanaiLab/eitan/UTRome/';
genesFiles = dir(strcat(workDir,'depth/*.txt'));
genesFiles = {genesFiles.name};
sampleDirs = dir(strcat(workDir,'samples/*'));
samples = {sampleDirs.name}; samples = samples(3:end);
cellAppearenceThreshold = 2;
readThreshold = 20;
minDistanceBetweenPopulations = 20;
outFile = fopen('isoformUsage_whole_embryos.csv','w');
fprintf(outFile,'%s,%s,','gene_name','peak_Location');
for n=1:length(samples)
     fprintf(outFile,'%s,',samples{n});
end
fprintf(outFile,'\n');

for geneIdx=1:size(genesFiles,2)
    
    gene = regexprep(genesFiles{geneIdx},'.txt','');
    [gene]
    file = strcat(workDir,'depth/',genesFiles{geneIdx});
    listing = dir(file);
    if (listing.bytes > 0)
        A = dlmread(file,'\t');
    
    if (sum(logical(sum(A)>readThreshold)) >= cellAppearenceThreshold) % at least 100 reads in 2 samples
            
        %generating a histogram based on all samples (will be skewed to highly represented
        smoothedCumulativeHistogram = tsmovavg(sum(A,2)','s', minDistanceBetweenPopulations);
        heightThreshold = max(5, sum(sum(A,2)) * 2e-3);
        [peakHeights,peakLocations] = findpeaks( smoothedCumulativeHistogram,'MINPEAKDISTANCE', minDistanceBetweenPopulations, ...
                                                'MINPEAKHEIGHT',heightThreshold);
        
        if (size(peakLocations,2) > 0)
            peakLocations = peakLocations - minDistanceBetweenPopulations / 2; %centralizing the peaks
            
            % for each peak
            for peakIdx=1:size(peakLocations,2)
                fprintf(outFile,'%s,%d,',gene,peakLocations(peakIdx));
                % for each sample
                for i=1:size(A,2)
                    sampleHistogram = A(:,i);
                    readNum = sum ( sampleHistogram ( ...
                         peakLocations(peakIdx) - minDistanceBetweenPopulations / 2 : ...
                         peakLocations(peakIdx) + minDistanceBetweenPopulations / 2 ) );
                    fprintf(outFile,'%d,',readNum);
                end
                fprintf(outFile,'\n');
            end
            
            % for all peaks
            fprintf(outFile,'%s,total,',gene);
            for i=1:size(A,2)
                 sampleHistogram = A(:,i);
                 readNum = sum(sampleHistogram);
                 fprintf(outFile,'%d,',readNum);
            end
            fprintf(outFile,'\n');
        end
    end
    end
end
fclose('all');

