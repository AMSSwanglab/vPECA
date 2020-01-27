function p_values=coeff_t_test(y,y_hat,X,beta)
n=length(y);
mean_X=mean(X);
p_values=ones(size(X,2),1);
for i=1:size(X,2)
    if beta(i)~=0
        SE=sqrt((y-y_hat)'*(y-y_hat)/(n-2))/sqrt((X(:,i)-mean(X(:,i)))'*(X(:,i)-mean(X(:,i))));
        t_n_2=beta(i)/SE;
        p_values(i)=2-2*tcdf(abs(t_n_2),n-2);
    end             
end
end