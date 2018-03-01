% Copied from gampmatlab toolbox
% http://gampmatlab.wikia.com/wiki/Generalized_Approximate_Message_Passing
%
% Finds the flipped, circ-shifted, and phase-rot version of "in" that best matches "ref"
% 
% [out,fxnvec,fxnmat]=disambig2Drfft(in,ref,n1,n2)

function [out,fxnvec,fxnmat]=disambig2Drfft(in,ref,n1,n2)

 % extract input size
 [n1in,n2in]=size(in);

 % convert to rectangular
 in=reshape(in,n1,n2);
 ref=reshape(ref,n1,n2);

 out = nan(size(in));
 minErr = inf;
 for flip=0:1
   if flip, 
     in_flip = fliplr(flipud(in)); 
   else
     in_flip = in;
   end
   for kk=0:n1-1
     for jj=0:n2-1
       in_shift = circshift(in_flip,[kk jj]);
       angl = sign(sum(sum(conj(in_shift).*ref)));
       in_rot = in_shift*angl;
       err = norm(in_rot-ref,'fro');
       if err<minErr
         minErr = err;
         out = in_rot; 
         flip_best = flip;
         kk_best = kk;
         jj_best = jj;
         angl_best = angl;
       end
     end
   end
 end

 % restore input size
 out = reshape(out,n1in,n2in);

 if nargout>1
   if flip_best
     fxnvec = @(invec) reshape(circshift(fliplr(flipud(reshape(invec,n1,n2))),[kk_best jj_best])*angl_best,n1in,n2in);
     fxnmat = @(inmat) squeeze(reshape(circshift(flipdim(flipdim(reshape(inmat,n1,n2,size(inmat,2)),1),2),[kk_best jj_best])*angl_best,n1in,n2in,size(inmat,2)));
   else
     fxnvec = @(invec) reshape(circshift(reshape(invec,n1,n2),[kk_best jj_best])*angl_best,n1in,n2in);
     fxnmat = @(inmat) squeeze(reshape(circshift(reshape(inmat,n1,n2,size(inmat,2)),[kk_best jj_best])*angl_best,n1in,n2in,size(inmat,2)));
   end
 end
