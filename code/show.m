function [] = show(img, new_img, recon_img, xx, yy, win)
% function [] = show(img, new_img, recon_img, video)
    maxo = max(img(:));
    mino = min(img(:));
    
%     outputVideo = VideoWriter(fullfile('../../Dropbox/videos',video), 'MPEG-4');
%     outputVideo.FrameRate = 3;
%     outputVideo.Quality = 70;
%     open(outputVideo)
% %     outputVideo1 = VideoWriter(fullfile('.','ps_output.avi'));
% %     open(outputVideo1)
    
    new_img1 = (new_img - mino)./(maxo - mino);
    recon_img1 = (recon_img - mino)./(maxo - mino);
    figure ('Position', [100 100 600 240]);
%     figure (1);
%     axis equal;
%     axis tight;
    RMSE = norm((recon_img(:) - img(:)),2)/norm(img(:),2);
    
    for i = 1:size(img,3)
    
        
        delete(findall(gcf,'type','annotation'))
        h(1) = subplot(1,3,1);
        set(h(1), 'position', [0.062, 0.225, 0.25, 0.7] );
        imagesc(img(:,:,i), [0, maxo]); colormap('jet'); 
        impixelinfo;
        title('Original');
        axis equal;
        axis tight;
        h(1) = subplot(1,3,2);
        set(h(1), 'position', [0.389, 0.225, 0.25, 0.7] );
        imagesc(new_img(:,:,i), [0, maxo]); colormap('jet'); 
        impixelinfo;
        axis equal;
        axis tight;
        title('Undersampled');
        h(1) = subplot(1,3,3);
        set(h(1), 'position', [0.712, 0.225, 0.25, 0.7] );
        imagesc(recon_img(:,:,i), [0, maxo]); colormap('jet'); 
        impixelinfo;
        title('Reconstructed');
        axis equal;
        axis tight;
%         title(strcat('Channel- ', num2str(i)));
        cb = colorbar('Location', 'southoutside');
        a = get(cb);
        a.Position;
        annotation('textbox', [0.05, 0.075, 0.15, 0.1], 'String', strcat('Channel- ', num2str(i)));
        set(cb,'Position',[0.238 0.1 0.488 0.05]);
        annotation('textbox', [0.763, 0.075, 0.2, 0.09], 'String', strcat('RMSE: ', num2str(RMSE)));
        
%         writeVideo(outputVideo,getframe(gcf));
% %         writeVideo(outputVideo1,recon_img1(:,:,i));
        pause(0.2);
    end
    
%     close(outputVideo);
% %     close(outputVideo1);
    figure;
    recon_spectra = zeros(1,size(img,3));
    orig_spectra = zeros(1,size(img,3));
    new_spectra = zeros(1,size(img,3));
    spectral_ind = zeros(1,size(img,3));
    
    for i=1:size(img,3)
        recon_pix = sum(sum(recon_img(xx:min(xx+win-1,size(img,1)),yy:min(yy+win-1,size(img,2)),i)))/(win*win);
        orig_pix = sum(sum(img(xx:min(xx+win-1,size(img,1)),yy:min(yy+win-1,size(img,2)),i)))/(win*win);
        new_pix = sum(sum(new_img(xx:min(xx+win-1,size(img,1)),yy:min(yy+win-1,size(img,2)),i)))/(win*win);
        recon_spectra(:,i) = recon_pix;
        orig_spectra(:,i) = orig_pix;
        new_spectra(:,i) = new_pix;
        spectral_ind(:,i) = i;
    end
        
    plot(spectral_ind, recon_spectra);
    hold on;
    plot(spectral_ind, orig_spectra);
    hold on;
    plot(spectral_ind, new_spectra);
    hold on;
    legend({'reconstructed','original'},'Location','northeast')
%     title(strcat('Spectral plot at pixel (', num2str(xx), ', ', num2str(yy), ') for',{' '}, video, ' with', {' '}, num2str(fh), '% sampling'), 'fontsize', 9);
%     title(strcat('Spectral plot at pixel (', num2str(xx), ', ', num2str(yy), ') for',{' '}, video, ' with every', {' '}, num2str(fh), 'th pixel sampling'), 'fontsize', 9);
        xlabel('Spectral band index');
        ylabel('Intensity values');
%     saveas(gcf, strcat(plotname, '.png'));
%         legend({'recon tissue', 'pureParaffin', 'mixed'},'Location','northwest')
%     end
% end
