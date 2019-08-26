% ad-hoc approach to spectral separation using peak identification and explicit scaling
function [height_derivatives, zero_crossings, boundaries, x] = spectral_id(spectra, pureParaffin)
   N_lambda = size(spectra, 3);
   figure(1);
   spectra(1,1,:) = spectra - min(spectra(1,1,:));
   
   [max_height, ~] = max(spectra(1,1,:));
   
   height_threshold = 0.25 * max_height;
   
   height_derivatives = zeros(size(spectra,3),1);
   height_derivatives(1,1) = 0;
   for i=2:N_lambda
       height_derivatives(i,1) = spectra(1,1,i) - spectra(1,1,i-1);
   end
   
   [max_height_derivative, ~] = max(height_derivatives(:,1));
   [min_height_derivative, ~] = min(height_derivatives(:,1));
   height_derivative_thresholdp = 0.05 * max_height_derivative;
   height_derivative_thresholdn = 0.05 * min_height_derivative;
   
   zero_crossings = zeros(N_lambda,1);
   for i=2:N_lambda-1
      if height_derivatives(i,1) > 0  && height_derivatives(i+1,1) < 0 && (spectra(1,1,i+1)+spectra(1,1,i))/2 > height_threshold
         zero_crossings(i,1) = 1;
      end
   end
   
   peaks = find(zero_crossings==1);
   peaks_iter = peaks';
   boundaries = zeros(numel(peaks_iter),2);
   ctr1 = 1;
   for i=peaks_iter
     ctr = 1;
     while height_derivatives(i+ctr,1) < height_derivative_thresholdn && height_derivatives(i-ctr,1) > height_derivative_thresholdp && i+ctr<=N_lambda && i-ctr>=1 
         ctr = ctr + 1;
     end
     boundaries(ctr1,1) = i-ctr;
     boundaries(ctr1,2) = i+ctr;
     ctr1 = ctr1 + 1;
   end
   
   indices = peaks;
   indices1 = boundaries(:);
   
%    scatter(indices1,squeeze(spectra(1,1,indices1)), '+');
%    hold on;
%    scatter(indices,squeeze(spectra(1,1,indices)), 'o');
%    hold on; 
%    plot(1:N_lambda, height_derivatives);
%    hold on;
   plot(1:N_lambda, squeeze(spectra(1,1,:)));
   hold on;
%    plot(1:N_lambda, zeros(N_lambda,1));
   
   cmax = zeros(numel(peaks_iter),1);
   atom_idx = ones(numel(peaks_iter),2);    
   peak_left_end = zeros(numel(peaks_iter),1);
   for i1=1:size(pureParaffin,1)
       for j1 = 1:size(pureParaffin,2)
            x = squeeze(pureParaffin(i1,j1,:));
            x = x - min(x);
%             plot(1:N_lambda, squeeze(x));
%             hold on;
            for k1 = 1:numel(peaks_iter)
                for shift=-5:5
                  if shift + boundaries(k1,2) > N_lambda || shift + boundaries(k1,1) < 1 
                      continue;
                  end
                  if (k1 < numel(peaks_iter) && shift+boundaries(k1,2) > boundaries(k1+1,1)) || (k1 > 1 && shift + boundaries(k1,1) < boundaries(k1-1,2))
                      continue;
                  end
                  c = NCC(squeeze(spectra(1,1,boundaries(k1,1):boundaries(k1,2))), x(boundaries(k1,1)+shift:boundaries(k1,2)+shift));
                  if (c > cmax(k1,1))
%                       disp("Updating...");
%                       figure;
%                       plot(1:size(squeeze(spectra(1,1,boundaries(k1,1):boundaries(k1,2))),1), squeeze(spectra(1,1,boundaries(k1,1):boundaries(k1,2))));
%                       hold on;
%                       plot(1:size(squeeze(spectra(1,1,boundaries(k1,1):boundaries(k1,2))),1), x(boundaries(k1,1)+shift:boundaries(k1,2)+shift));
%                       hold on;    
                      cmax(k1,1) = c;
                      atom_idx(k1,:) = [i1 j1];
                      peak_left_end(k1,1) = shift + boundaries(k1,1);
                  end
                end
            end
%        plot(1:N_lambda, squeeze(x));
%        hold on;
       end
   end
   
   for k1=1:numel(peaks_iter)
        peak_left_end(k1);
        boundaries(k1,:);
        bound = min(peak_left_end(k1)+boundaries(k1,2)-boundaries(k1,1), N_lambda);
        spectra(1,1,peak_left_end(k1):bound) = spectra(1,1,boundaries(k1,1):boundaries(k1,1)+bound - peak_left_end(k1));
        if (peak_left_end(k1) - boundaries(k1,1) > 0)
           spectra(1,1,boundaries(k1,1):peak_left_end(k1)-1) = spline([boundaries(k1,1)-1 peak_left_end(k1)],[spectra(1,1,boundaries(k1,1)-1) spectra(1,1,peak_left_end(k1))],boundaries(k1,1):peak_left_end(k1)-1); 
        elseif (peak_left_end(k1) - boundaries(k1,1) < 0)
            spectra(1,1,peak_left_end(k1)+boundaries(k1,2)-boundaries(k1,1)+1:boundaries(k1,2)) = spline([peak_left_end(k1)+boundaries(k1,2)-boundaries(k1,1) boundaries(k1,2)+1],[spectra(1,1,peak_left_end(k1)+boundaries(k1,2)-boundaries(k1,1)) spectra(1,1,boundaries(k1,2)+1)],peak_left_end(k1)+boundaries(k1,2)-boundaries(k1,1)+1:boundaries(k1,2));
        end
%         plot(peak_left_end(k1):peak_left_end(k1)+boundaries(k1,2)-boundaries(k1,1), squeeze(pureParaffin(atom_idx(k1,1), atom_idx(k1,2), peak_left_end(k1):peak_left_end(k1)+boundaries(k1,2)-boundaries(k1,1)) - min(pureParaffin(atom_idx(k1,1), atom_idx(k1,2), peak_left_end(k1):peak_left_end(k1)+boundaries(k1,2)-boundaries(k1,1)))));
%         hold on;
   end
   
   plot(1:N_lambda, squeeze(spectra(1,1,:)));
   hold on;

%     scatter(peak_left_end(k1,1), 1e6 * squeeze(A(1,1,peak_left_end(k1,1), atom_idx(k1,1))), '+');
%     line([peak_left_end peak_left_end], get(gca, 'ylim'));
    
%    y = sgolayfilt(x,20,207); % No need of filtered baseline removal
%    plot(1:N_lambda, squeeze(y));

    p = squeeze(spectra(1,1,:));
    pval = zeros(size(squeeze(spectra(1,1,:)),1), numel(peaks_iter));
    for i2=1:numel(peaks_iter)
        pval(:,i2) = squeeze(pureParaffin(atom_idx(i2,1), atom_idx(i2,2),:));
    end
    
    coeff_matrix = zeros(numel(peaks_iter));
    for i2=1:size(coeff_matrix,1)
        for j2=1:size(coeff_matrix,2)
            coeff_matrix(i2,j2) = sum(pval(:,i2).*pval(:,j2));
        end
    end
    
    rhs = zeros(numel(peaks_iter),1);
    for i2=1:numel(peaks_iter)
        rhs(i2,1) = sum(p.*pval(:,i2));
    end
    
    scalers = coeff_matrix\rhs;
    
    tmp_spec = squeeze(spectra(1,1,:));
    for i2=1:numel(peaks_iter)
        tmp_spec = tmp_spec - scalers(i2).*pval(:,i2);
    end
    spectra(1,1,:) = reshape(tmp_spec, [1 N_lambda]);
    
    spectra(spectra(1,1,:)<0) = 0;
    plot(1:N_lambda, squeeze(spectra(1,1,:)));
    hold on;
    
    legend({'raw spectra', 'shifted spectra', 'subtracted spectra'},'Location','northeast')

end

function val = NCC(vx, vy) % Here, vx is the template
    if size(vx,2) > 1 || size(vy,2) > 1
        disp("NOT 1D vectors");
    end
    prod = (vx - mean(vx)) .* (vy - mean(vy));
    val = sum(prod)/(norm(vx - mean(vx),2)*norm(vy - mean(vy),2));
    n = numel(vx); % Considering vx is 1D
    val = val/n;
end