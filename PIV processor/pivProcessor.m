%% Made by Nout van den Bos
%This function requires you to have extrema.m and extrema2.m in this
%directory. These are opensource functions publicly available.

function pivProcessor(foldread,filestore,fileread,ws,ovlap,dt,pix_size,...
                      xlength,ylength,M,x0,y0,SNRmin,filterSignal)


    %read the image and divide into 2 images
    if foldread
        im = imread([foldread fileread]);
    else
        im = imread(fileread);
    end

    pic1 = im(1:ylength,:);
    pic2 = im(ylength+1:end,:);

    %determine the number of windows you can fit in the pictures
    xwindows = floor(xlength/ws);
    ywindows = floor(ylength/ws);
    
    %declare the lists of the output variables
    u = zeros(1,xwindows*ywindows);
    v = zeros(1,xwindows*ywindows);%[];
    xpix = zeros(1,xwindows*ywindows);% [];
    ypix = zeros(1,xwindows*ywindows);% [];
    SNR = zeros(1,xwindows*ywindows);% [];
    
    %this is the index for the results
    index = 0;
    
    for i = 1:ywindows
        for j = 1:xwindows 
            
            index = index+1;
            %calculate the indexes for the pictures, take into account the
            %overlap percentage.
            idxstartY = max((i-1)*ws+1-round(ws*ovlap),1);
            idxendY   = idxstartY+ws-1;
            idxstartX = max((j-1)*ws+1-round(ws*ovlap),1);
            idxendX   = idxstartX+ws-1;
        
            %establish the windows, convert to double and substract the
            %mean in order to get better separated peaks in the cross
            %correlation function.
            wind1 = double(pic1(idxstartY:idxendY,idxstartX:idxendX));
            wind1 = wind1- mean(wind1,'all');  
            wind2 = double(pic2(idxstartY:idxendY,idxstartX:idxendX));
            wind2 = wind2-mean(wind2,'all');

            %calculate the correlations, find the position of the maximum
            c = xcorr2(wind1,wind2);
            [ypeak,xpeak]    = find(c==max(c(:)));
     
            %locate values of all peaks, in order to calculate the signal
            %to noise ratio
            if (filterSignal)
                [allpeaks,~] = extrema2(c);
                SNRtemp = allpeaks(1)/allpeaks(2);          
            else
                SNRtemp = 1;
            end
            
            %Use gaussian fit for sub-pixel precision
            [ypeak, xpeak] = gaussianfit(c,ypeak,xpeak);

            %calculate offsets and velocities etc
            yoffSet = (ws-ypeak);
            xoffSet = (ws-xpeak);
            
            %attach the results to the results vectors
            u(index) = xoffSet*pix_size/(M*dt);
            v(index) = -yoffSet*pix_size/(M*dt);
            xpix(index) = idxstartX;
            ypix(index) = idxstartY;
            SNR(index)  = SNRtemp;
            
        end 
    end

    %map pixel positions onto real positions
    x = (xpix-x0).*pix_size./M.*10^-3;
    y = -(ypix-y0).*pix_size./M.*10^-3;
    
    %Only use masks if user asks it
    if(filterSignal)
        
    %Signal to Noise mask
        mask1 = (SNR<SNRmin) ;

        %outliers mask
        mask2 = isoutlier(u,'movmedian',10);
        mask3 = isoutlier(v,'movmedian',10);

        %total mask
        mask = mask1|mask2|mask3;

        %inverse it as to follow the conventions
        mask = ~mask;
        
    else
        mask = ones(1,length(SNR));
    end
    
    %write in .dat format
    results = [x; y; u; v; mask]; 
    fileID = fopen(filestore, 'w');
    fprintf(fileID,'%6s %6s %6s %6s %6s \n','x','y','u','v','mask');
    fprintf(fileID,'%4.4f %4.4f %2.6f %2.6f %1d\n',results );
    fclose(fileID);
end
