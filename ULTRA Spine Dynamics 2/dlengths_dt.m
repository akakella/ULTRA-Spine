function dlengths_dt = dlengths_dt(x2,y2,z2,T2,G2,P2,dx2,dy2,dz2,dT2,dG2,dP2)
%DLENGTHS_DT
%    DLENGTHS_DT = DLENGTHS_DT(X2,Y2,Z2,T2,G2,P2,DX2,DY2,DZ2,DT2,DG2,DP2)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    10-Dec-2015 16:13:58

t2 = sqrt(3.0);
t3 = sin(G2);
t4 = cos(G2);
t5 = cos(T2);
t12 = z2.*4.0e1;
t13 = t2.*t3.*3.0;
t14 = t4.*t5.*3.0;
t6 = t12-t13-t14+3.0;
t7 = cos(P2);
t8 = sin(P2);
t9 = sin(T2);
t19 = x2.*4.0e1;
t20 = t2.*3.0;
t21 = t8.*t9.*3.0;
t22 = t2.*t4.*t7.*3.0;
t23 = t3.*t5.*t7.*3.0;
t10 = -t19+t20+t21-t22+t23;
t26 = y2.*4.0e1;
t27 = t7.*t9.*3.0;
t28 = t2.*t4.*t8.*3.0;
t29 = t3.*t5.*t8.*3.0;
t11 = t26-t27-t28+t29;
t15 = t12+t13-t14+3.0;
t16 = t2.*(3.0./4.0e1);
t17 = t8.*t9.*(3.0./4.0e1);
t18 = t3.*t5.*t7.*(3.0./4.0e1);
t24 = -t19-t20+t21+t22+t23;
t25 = t3.*t5.*t8.*(3.0./4.0e1);
t30 = t26-t27+t28+t29;
t31 = dz2.*1.2e2;
t32 = dx2.*x2.*1.6e3;
t33 = dy2.*y2.*1.6e3;
t34 = dz2.*z2.*1.6e3;
t35 = dG2.*t3.*t7.*2.7e1;
t36 = dP2.*t4.*t8.*2.7e1;
t37 = dG2.*t3.*t5.*9.0;
t38 = dT2.*t4.*t9.*9.0;
t39 = dx2.*t2.*t4.*t7.*1.2e2;
t40 = dP2.*t2.*t7.*t9.*9.0;
t41 = dT2.*t2.*t5.*t8.*9.0;
t42 = dG2.*t3.*t5.*z2.*1.2e2;
t43 = dT2.*t4.*t9.*z2.*1.2e2;
t44 = dP2.*t8.*t9.*y2.*1.2e2;
t45 = dy2.*t3.*t5.*t8.*1.2e2;
t46 = dG2.*t2.*t3.*t8.*y2.*1.2e2;
t47 = dG2.*t2.*t4.*t5.*t7.*9.0;
t48 = dG2.*t4.*t5.*t8.*y2.*1.2e2;
t49 = dP2.*t3.*t5.*t7.*y2.*1.2e2;
t50 = dP2.*t3.*t5.*t8.*x2.*1.2e2;
t51 = dT2.*t3.*t7.*t9.*x2.*1.2e2;
t62 = t2.*t5.*t8.*(3.0./4.0e1);
t63 = t2.*t3.*t7.*t9.*(3.0./4.0e1);
t64 = t17+t18+t62-t63+x2;
t52 = abs(t64);
t59 = t4.*t5.*(3.0./4.0e1);
t60 = t2.*t4.*t9.*(3.0./4.0e1);
t61 = t59-t60+z2-3.0./4.0e1;
t53 = abs(t61);
t55 = t7.*t9.*(3.0./4.0e1);
t56 = t2.*t5.*t7.*(3.0./4.0e1);
t57 = t2.*t3.*t8.*t9.*(3.0./4.0e1);
t58 = -t16-t25+t55+t56+t57+y2;
t54 = abs(t58);
t72 = t17+t18-t62+t63+x2;
t65 = abs(t72);
t69 = t59+t60+z2-3.0./4.0e1;
t66 = abs(t69);
t68 = t16-t25+t55-t56-t57+y2;
t67 = abs(t68);
t70 = dP2.*t5.*t8.*2.7e1;
t71 = dT2.*t7.*t9.*2.7e1;
t73 = dT2.*t2.*t4.*t5.*9.0;
t74 = dT2.*t5.*t7.*y2.*1.2e2;
t75 = dP2.*t2.*t8.*t9.*9.0;
t76 = dP2.*t7.*t9.*x2.*1.2e2;
t77 = dT2.*t5.*t8.*x2.*1.2e2;
t78 = dP2.*t2.*t5.*t7.*x2.*1.2e2;
t79 = dG2.*t2.*t3.*t9.*z2.*1.2e2;
t80 = dG2.*t2.*t4.*t5.*t8.*9.0;
t81 = dP2.*t2.*t3.*t5.*t7.*9.0;
t82 = dG2.*t4.*t5.*t7.*x2.*1.2e2;
t83 = dT2.*t3.*t8.*t9.*y2.*1.2e2;
t84 = dG2.*t2.*t4.*t8.*t9.*y2.*1.2e2;
t85 = dP2.*t2.*t3.*t7.*t9.*y2.*1.2e2;
t86 = dT2.*t2.*t3.*t5.*t8.*y2.*1.2e2;
t87 = dP2.*t2.*t3.*t8.*t9.*x2.*1.2e2;
t88 = t2.*t3.*(3.0./4.0e1);
t89 = -t12+t13+t14+3.0;
t90 = t2.*t4.*t8.*(3.0./4.0e1);
t91 = t20-t26+t27+t28-t29;
t92 = t2.*t4.*t7.*(3.0./4.0e1);
t93 = t19-t21+t22-t23;
t94 = dG2.*t2.*t4.*9.0;
t95 = dz2.*t2.*t3.*1.2e2;
t96 = dG2.*t2.*t4.*z2.*1.2e2;
t97 = dT2.*t2.*t5.*t7.*9.0;
t98 = dy2.*t2.*t4.*t8.*1.2e2;
t99 = dG2.*t2.*t3.*t7.*x2.*1.2e2;
t100 = dP2.*t2.*t4.*t8.*x2.*1.2e2;
t101 = dT2.*t2.*t3.*t8.*t9.*9.0;
t102 = dP2.*t2.*t4.*t7.*y2.*1.2e2;
t103 = t20+t26-t27-t28+t29;
t104 = t59+t88-z2+3.0./4.0e1;
t105 = sign(t104);
t106 = 1.0./t105.^2;
t107 = t89.^2;
t108 = t106.*t107.*6.25e-4;
t109 = -t17-t18+t92+x2;
t110 = sign(t109);
t111 = 1.0./t110.^2;
t112 = t93.^2;
t113 = t111.*t112.*6.25e-4;
t114 = dy2.*t2.*1.2e2;
t115 = dz2.*t4.*t5.*1.2e2;
t116 = dy2.*t7.*t9.*1.2e2;
t117 = dG2.*t3.*t8.*2.7e1;
t118 = dx2.*t8.*t9.*1.2e2;
t119 = dx2.*t3.*t5.*t7.*1.2e2;
t120 = -t19+t21+t22+t23;
t121 = t12+t13-t14-3.0;
t122 = -t20+t26-t27+t28+t29;
t123 = dP2.*t4.*t7.*2.7e1;
t124 = t17+t18+t92-x2;
t125 = sign(t124);
t126 = 1.0./t125.^2;
t127 = t120.^2;
t128 = t126.*t127.*6.25e-4;
t129 = -t59+t88+z2-3.0./4.0e1;
t130 = sign(t129);
t131 = 1.0./t130.^2;
t132 = t121.^2;
t133 = t131.*t132.*6.25e-4;
t134 = t20+t26-t27+t28+t29;
dlengths_dt = [1.0./sqrt(t10.^2.*1.0./sign(t16+t17+t18-x2-t2.*t4.*t7.*(3.0./4.0e1)).^2.*6.25e-4+t11.^2.*1.0./sign(t25+y2-t7.*t9.*(3.0./4.0e1)-t2.*t4.*t8.*(3.0./4.0e1)).^2.*6.25e-4+t6.^2.*1.0./sign(z2-t2.*t3.*(3.0./4.0e1)-t4.*t5.*(3.0./4.0e1)+3.0./4.0e1).^2.*6.25e-4).*(t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51-dx2.*t2.*1.2e2-dG2.*t2.*t4.*9.0-dx2.*t8.*t9.*1.2e2-dy2.*t7.*t9.*1.2e2-dz2.*t2.*t3.*1.2e2-dz2.*t4.*t5.*1.2e2-dx2.*t3.*t5.*t7.*1.2e2-dy2.*t2.*t4.*t8.*1.2e2-dP2.*t7.*t9.*x2.*1.2e2-dT2.*t5.*t8.*x2.*1.2e2-dT2.*t5.*t7.*y2.*1.2e2-dG2.*t2.*t4.*z2.*1.2e2-dP2.*t2.*t3.*t5.*t8.*9.0-dT2.*t2.*t3.*t7.*t9.*9.0-dG2.*t2.*t3.*t7.*x2.*1.2e2-dG2.*t4.*t5.*t7.*x2.*1.2e2-dP2.*t2.*t4.*t8.*x2.*1.2e2-dP2.*t2.*t4.*t7.*y2.*1.2e2-dT2.*t3.*t8.*t9.*y2.*1.2e2).*6.25e-4,1.0./sqrt(t30.^2.*1.0./sign(t25+t90+y2-t7.*t9.*(3.0./4.0e1)).^2.*6.25e-4+t24.^2.*1.0./sign(-t16+t17+t18+t92-x2).^2.*6.25e-4+t15.^2.*1.0./sign(t88+z2-t4.*t5.*(3.0./4.0e1)+3.0./4.0e1).^2.*6.25e-4).*(t31+t32+t33+t34+t35+t36+t37+t38-t39-t40-t41+t42+t43+t44+t45-t46-t47+t48+t49+t50+t51+t94+t95+t96+t98+t99+t100+t102+dx2.*t2.*1.2e2-dx2.*t8.*t9.*1.2e2-dy2.*t7.*t9.*1.2e2-dz2.*t4.*t5.*1.2e2-dx2.*t3.*t5.*t7.*1.2e2-dP2.*t7.*t9.*x2.*1.2e2-dT2.*t5.*t8.*x2.*1.2e2-dT2.*t5.*t7.*y2.*1.2e2+dP2.*t2.*t3.*t5.*t8.*9.0+dT2.*t2.*t3.*t7.*t9.*9.0-dG2.*t4.*t5.*t7.*x2.*1.2e2-dT2.*t3.*t8.*t9.*y2.*1.2e2).*6.25e-4,1.0./sqrt(t52.^2+t53.^2+t54.^2).*(t37+t38-t42-t43-t44-t48-t49-t50-t51+t70+t71+t73+t74+t75+t76+t77+t78+t79+t80+t81+t82+t83+t84+t85+t86+t87+dx2.*t52.*sign(t64).*1.6e3+dy2.*t54.*sign(t58).*1.6e3+dz2.*t53.*sign(t61).*1.6e3-dG2.*t2.*t3.*t9.*9.0-dG2.*t4.*t8.*t9.*2.7e1-dP2.*t3.*t7.*t9.*2.7e1-dT2.*t2.*t5.*t7.*9.0-dT2.*t3.*t5.*t8.*2.7e1-dT2.*t2.*t3.*t8.*t9.*9.0-dT2.*t2.*t8.*t9.*x2.*1.2e2-dP2.*t2.*t5.*t8.*y2.*1.2e2-dT2.*t2.*t7.*t9.*y2.*1.2e2-dT2.*t2.*t4.*t5.*z2.*1.2e2-dG2.*t2.*t4.*t7.*t9.*x2.*1.2e2-dT2.*t2.*t3.*t5.*t7.*x2.*1.2e2).*6.25e-4,1.0./sqrt(t65.^2+t66.^2+t67.^2).*(t37+t38-t42-t43-t44-t48-t49-t50-t51+t70+t71-t73+t74-t75+t76+t77-t78-t79-t80-t81+t82+t83-t84-t85-t86-t87+t97+t101+dx2.*t65.*sign(t72).*1.6e3+dy2.*t67.*sign(t68).*1.6e3+dz2.*t66.*sign(t69).*1.6e3+dG2.*t2.*t3.*t9.*9.0-dG2.*t4.*t8.*t9.*2.7e1-dP2.*t3.*t7.*t9.*2.7e1-dT2.*t3.*t5.*t8.*2.7e1+dT2.*t2.*t8.*t9.*x2.*1.2e2+dP2.*t2.*t5.*t8.*y2.*1.2e2+dT2.*t2.*t7.*t9.*y2.*1.2e2+dT2.*t2.*t4.*t5.*z2.*1.2e2+dG2.*t2.*t4.*t7.*t9.*x2.*1.2e2+dT2.*t2.*t3.*t5.*t7.*x2.*1.2e2).*6.25e-4,1.0./sqrt(t108+t113+t91.^2.*1.0./sign(t16-t25+t55+t90-y2).^2.*6.25e-4).*(t31-t32-t33-t34+t37+t38-t39-t42-t43-t44-t45-t46-t48-t49-t50-t51+t74+t75+t76+t77+t80+t81+t82+t83-t94+t95+t96-t97+t98+t99+t100-t101+t102+t114+t115+t116+t117+t118+t119-dP2.*t4.*t7.*2.7e1).*(-6.25e-4),1.0./sqrt(t108+t113+t103.^2.*1.0./sign(t16+t25-t55-t90+y2).^2.*6.25e-4).*(t31-t32-t33-t34+t37+t38-t39-t42-t43-t44-t45-t46-t48-t49-t50-t51+t74-t75+t76+t77-t80-t81+t82+t83-t94+t95+t96+t97+t98+t99+t100+t101+t102-t114+t115+t116-t117+t118+t119+t123).*(-6.25e-4),1.0./sqrt(t128+t133+t122.^2.*1.0./sign(-t16+t25-t55+t90+y2).^2.*6.25e-4).*(t31-t32-t33-t34+t37+t38+t39-t42-t43-t44-t45+t46-t48-t49-t50-t51+t74+t75+t76+t77+t80+t81+t82+t83+t94-t95-t96-t97-t98-t99-t100-t101-t102+t114+t115+t116-t117+t118+t119+t123).*(-6.25e-4),1.0./sqrt(t128+t133+t134.^2.*1.0./sign(t16+t25-t55+t90+y2).^2.*6.25e-4).*(-t31+t32+t33+t34-t37-t38-t39+t42+t43+t44+t45-t46+t48+t49+t50+t51-t74+t75-t76-t77+t80+t81-t82-t83-t94+t95+t96-t97+t98+t99+t100-t101+t102+t114-t115-t116-t117-t118-t119+t123).*6.25e-4];
