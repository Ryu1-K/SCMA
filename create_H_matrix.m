function H = create_H_matrix( CH, SC )
    % Create fading matrix
    tmp_fad = ( randn(SC.F, SC.N*SC.L) + 1i*randn(SC.F, SC.N*SC.L) )*sqrt(CH.Nh/2);
%     tmp_fad = ones(SC.F, SC.N*SC.L);
    
    CH.fad =  spalloc(SC.F*SC.N, SC.F*SC.L, SC.F*SC.N*SC.L);
    for l = 1:SC.L
        for n =1:SC.N
            CH.fad((n-1)*SC.F+1 : n*SC.F ,(l-1)*SC.F+1 : l*SC.F) = diag( tmp_fad(:,l+(n-1)*SC.F) );
        end
    end

    H = CH.fad*SC.V*SC.A;
%     H = CH.fad*SC.VM_matrix{randi(SC.C)};
end

