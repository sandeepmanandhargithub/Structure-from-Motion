function imagePoints = world2image(T, X)
n = size(X);

imagePoints = [];
    for i =1:n
        imagePoints = [imagePoints T*X(i,:)'];
    end
     for i =1:length(imagePoints)
        imagePoints(:,i) = imagePoints(:,i)/imagePoints(3,i);
    end
    
    
