function [a, Ihat, error] = sparsify_bruno(Phi, I)

    [N, M, w] = size(Phi);
    T = size(I, 2);
    num_iter = 50;
    beta = 1;
    sigma = .5;
    bos = beta/sigma;
    eta = 1e-2;
    gain = sqrt(sum(sum(Phi.*Phi),3))';
    
    lambda_N = 1;
    marg = floor(w/2);
    num_valid = T - 2 * marg;
    sbuf = 1:marg;
    ebuf = T - marg + 1:T;
    
    I(:,sbuf)=0;
    I(:,ebuf)=0;
    
    error_record = zeros(N, num_iter);
    
    a = 0;
    for t = 1:w
        tt = t:T - (w - t);
        a = a + Phi(:, :, t)' * I(:, tt);
    end
    
    for i=1:M
        a(i,:) = a(i,:)/(gain(i)*gain(i));
    end
    
    Ihat = zeros(N, T);
    for t = 1:w
        tt = t:T - (w - t);
        Ihat(:, tt) = Ihat(:, tt) + Phi(:, :, t) * a;
    end
    
    error = I - Ihat;
    error(:,sbuf) = 0;
    error(:,ebuf) = 0;
    
    E=(0.5*lambda_N*sum(error(:).^2)+beta*sum(S(a(:)/sigma)))/num_valid;
    
    for i = 1:num_iter
        i
        da = zeros(size(a));
        for t = 1:w
            tt = t:T - (w - t);
            da = da + Phi(:, :, t)' * error(:, tt);
        end
        
        bospa = bos*Sp(a/sigma);
        %da = da - bospa;
        da = da/norm(da) - bos*Sp(a/sigma);
        a = a + eta * da;
        
        for t = 1:w
            tt = t:T - (w - t);
            Ihat(:, tt) = Ihat(:, tt) + Phi(:, :, t) * da;
        end

        error = I - Ihat;
        error(:,sbuf) = 0;
        error(:,ebuf) = 0;
        error_record(:, i) = sum(error.^2, 2)/num_valid;
        
        Eold=E;
        E= (0.5*lambda_N*sum(error(:).^2)+beta*sum(S(a(:)/sigma)))/num_valid;
        fprintf('\rE=%10.4f, dEpct=%6.2f%%',E,100*(E-Eold)/Eold);
%         if (E-Eold)>0
%             a=[];
%             return
%         end
        
    end
    
    figure(200)
    for i=1:16;
        subplot(4,4,i)
        plot(error_record(i,:))
    end
    suptitle('Coefficient Reconstruction Error');

    figure(300)
    for j = 1:16;
        subplot(4, 4, j)
        plot(a(j, :));
    end
    suptitle('Sparse Coefficients');
    end