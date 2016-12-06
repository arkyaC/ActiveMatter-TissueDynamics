function saveFrames(frames)
    skipFrames = 3;
    for i = 1:length(frames)
        if mod(i,skipFrames)~=0
            continue
        end
        thisframe=frames(1,i);
        thisImage=frame2im(thisframe);
        thisFile = sprintf('/home/jarvis/Documents/VertexModel/movie(noise=0.6,rho=0.7)/frames/frame_%04d.jpg',i/skipFrames);
        imwrite(thisImage,thisFile);
    end
end