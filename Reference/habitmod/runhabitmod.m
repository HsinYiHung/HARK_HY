dirnam='./'; % directory of this example 
parnam='habitmodparams';
modnam='habitmod';
[cof, scof, cofb, param_, eqname_,...
 endog_, eqtype_, vtype_, neq, nlag, nlead, rts, lgrts,aimcode]=...
SPSolve(dirnam, modnam, parnam);

