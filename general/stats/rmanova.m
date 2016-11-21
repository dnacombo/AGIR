function stats = rmanova(data,cov)

% stats = rmanova(data[,cov])

error('still needs work')

vectorize = @(x) x(:);

c = doconds(size(data));
Y = data(:);
X = cellfun(vectorize,c);

B = mvregress(X,Y);


