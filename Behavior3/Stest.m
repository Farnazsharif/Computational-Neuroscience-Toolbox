function [Pa,Pk] = Stest( x,y )


name1=repmat({'x'},1,length(x));
name2=repmat({'y'},1,length(y));

groups=[x,y];
names=[name1,name2];
[Pa,~,~] = anova1(groups, names,'off');
[Pk,~,~] = kruskalwallis(groups, names,'off');
% [~,Pt] = ttest(x,y);

end

% anova1:         tests the hypothesis that the samples in y are drawn from populations with the same mean 
%                 against the alternative hypothesis that the population means are not all the same.
%                 
%kruskalwallis:   returns the p-value for the null hypothesis that the data in each column of the matrix x comes from 
%                 the same distribution. The alternative hypothesis is that not all samples come from the same distribution.

%ttest:           returns a test decision for the null hypothesis that the data in x comes from a normal distribution
%                 with mean equal to zero and unknown variance.The alternative hypothesis is that the population distribution
%                 does not have a mean equal to zero