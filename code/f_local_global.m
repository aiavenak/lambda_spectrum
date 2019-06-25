function [G,D,B,steps] = f_local_global(SC,coor,global_local_params)

switch global_local_params
    case 'SPL_W_log'
        %fprintf(' SPL with Log transform \n');
        D = -log(SC);
        [G,steps,B] = f_FastFloyd_und(D);
        D(isinf(D)) = 0;
        
    case 'SPL_W_inv' 
        %fprintf(' SPL with Inv transform \n');
        D = 1./SC;
        [G,steps,B] = f_FastFloyd_und(D);
        D(isinf(D)) = 0;
   
    case 'bin' 
        %fprintf(' SPL with Inv transform \n');
        D = 1./(SC>0);
        [G,steps,B] = f_FastFloyd_und(D);
        D(isinf(D)) = 0;        
        
    case 'SPL_W_SI'
        %fprintf(' SPL based on Search Information \n');
        deg = diag(1./sum(SC,2));
        D = -log(deg*SC);
        [G,steps,B] = f_FastFloyd_und(D);
        D(isinf(D)) = 0;

    case 'ED_W_log'
        %fprintf(' SPL with Log transform \n');
        D = -log(SC);
        [~,steps,B] = f_FastFloyd_und(D);
        D(isinf(D)) = 0;
        G = squareform(pdist(coor));
        
    case 'ED_W_inv'
        %fprintf(' SPL with Inv transform \n');
        D = 1./SC;
        D(isinf(D)) = 0;
        G = squareform(pdist(coor));
        
    case 'ED_W_SI'
        %fprintf(' SPL based on Search Information \n');
        deg = diag(1./sum(SC,2));
        D = -log(deg*SC);
        D(isinf(D)) = 0;        
        G = squareform(pdist(coor));
        
    case 'SPL_W_inv_minus1'
        %fprintf(' Matching Index \n');
        %MI = matching_ind_und(SC);
        D = (1./SC)-1;
        [G,steps,B] = f_FastFloyd_und(D);
        D(isinf(D)) = 0;
        
        
end
