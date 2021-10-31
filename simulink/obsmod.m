function [ret,x0,str,ts,xts]=obsmod(t,x,u,flag);
%OBSMOD	is the M-file description of the SIMULINK system named OBSMOD.
%	The block-diagram can be displayed by typing: OBSMOD.
%
%	SYS=OBSMOD(T,X,U,FLAG) returns depending on FLAG certain
%	system values given time point, T, current state vector, X,
%	and input vector, U.
%	FLAG is used to indicate the type of output to be returned in SYS.
%
%	Setting FLAG=1 causes OBSMOD to return state derivatives, FLAG=2
%	discrete states, FLAG=3 system outputs and FLAG=4 next sample
%	time. For more information and other options see SFUNC.
%
%	Calling OBSMOD with a FLAG of zero:
%	[SIZES]=OBSMOD([],[],[],0),  returns a vector, SIZES, which
%	contains the sizes of the state vector and other parameters.
%		SIZES(1) number of states
%		SIZES(2) number of discrete states
%		SIZES(3) number of outputs
%		SIZES(4) number of inputs
%		SIZES(5) number of roots (currently unsupported)
%		SIZES(6) direct feedthrough flag
%		SIZES(7) number of sample times
%
%	For the definition of other parameters in SIZES, see SFUNC.
%	See also, TRIM, LINMOD, LINSIM, EULER, RK23, RK45, ADAMS, GEAR.

% Note: This M-file is only used for saving graphical information;
%       after the model is loaded into memory an internal model
%       representation is used.

% the system will take on the name of this mfile:
sys = mfilename;
new_system(sys)
simver(1.3)
if (0 == (nargin + nargout))
     set_param(sys,'Location',[135,253,1068,654])
     open_system(sys)
end;
set_param(sys,'algorithm',     'RK-45')
set_param(sys,'Start time',    '0.0')
set_param(sys,'Stop time',     '140')
set_param(sys,'Min step size', '0.01')
set_param(sys,'Max step size', '0.25')
set_param(sys,'Relative error','1e-3')
set_param(sys,'Return vars',   '')

add_block('built-in/Sum',[sys,'/','Sum6'])
set_param([sys,'/','Sum6'],...
		'orientation',3,...
		'inputs','-+',...
		'position',[361,315,394,330])

add_block('built-in/State-Space',[sys,'/','State-Space1'])
set_param([sys,'/','State-Space1'],...
		'A','Am',...
		'B','Bm',...
		'C','eye(n)',...
		'D','zeros(n,m)',...
		'position',[245,325,325,365])

add_block('built-in/To Workspace',[sys,'/','To Workspace8'])
set_param([sys,'/','To Workspace8'],...
		'mat-name','x',...
		'buffer','1000000',...
		'position',[485,65,535,85])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain3']])
set_param([sys,'/',['Matrix',13,'Gain3']],...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:',...
		'Mask Translate','K = @1;')
set_param([sys,'/',['Matrix',13,'Gain3']],...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','B\/',...
		'position',[580,45,610,75])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain1']])
set_param([sys,'/',['Matrix',13,'Gain1']],...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:',...
		'Mask Translate','K = @1;')
set_param([sys,'/',['Matrix',13,'Gain1']],...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','C\/',...
		'position',[400,85,430,115])

add_block('built-in/State-Space',[sys,'/','State-Space'])
set_param([sys,'/','State-Space'],...
		'A','A',...
		'B','B',...
		'C','eye(n)',...
		'D','zeros(n,m)',...
		'X0','x0',...
		'position',[280,79,355,121])

add_block('built-in/To Workspace',[sys,'/','To Workspace2'])
set_param([sys,'/','To Workspace2'],...
		'mat-name','xm',...
		'buffer','1000000',...
		'position',[430,362,480,378])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain12']])
set_param([sys,'/',['Matrix',13,'Gain12']],...
		'orientation',3,...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:')
set_param([sys,'/',['Matrix',13,'Gain12']],...
		'Mask Translate','K = @1;',...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','Lx\/',...
		'position',[400,205,430,235])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain11']])
set_param([sys,'/',['Matrix',13,'Gain11']],...
		'orientation',3,...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:')
set_param([sys,'/',['Matrix',13,'Gain11']],...
		'Mask Translate','K = @1;',...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','P*S\/',...
		'position',[325,250,355,280])

add_block('built-in/MATLAB Fcn',[sys,'/','Unv1'])
set_param([sys,'/','Unv1'],...
		'orientation',3,...
		'MATLAB Fcn','u/(norm(u)+delta)',...
		'Mask Display','UnV',...
		'Mask Type','unit vector',...
		'Mask Dialogue','Unit Vector|delta',...
		'Mask Translate','delta=@1;')
set_param([sys,'/','Unv1'],...
		'Mask Entries','0.001\/',...
		'position',[320,200,360,230])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain13']])
set_param([sys,'/',['Matrix',13,'Gain13']],...
		'orientation',2,...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:')
set_param([sys,'/',['Matrix',13,'Gain13']],...
		'Mask Translate','K = @1;',...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','L\/',...
		'position',[435,135,465,165])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain4']])
set_param([sys,'/',['Matrix',13,'Gain4']],...
		'orientation',2,...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:')
set_param([sys,'/',['Matrix',13,'Gain4']],...
		'Mask Translate','K = @1;',...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','C\/',...
		'position',[660,285,690,315])

add_block('built-in/To Workspace',[sys,'/','To Workspace3'])
set_param([sys,'/','To Workspace3'],...
		'orientation',3,...
		'mat-name','yobs',...
		'buffer','1000000',...
		'position',[520,265,570,285])

add_block('built-in/To Workspace',[sys,'/','To Workspace5'])
set_param([sys,'/','To Workspace5'],...
		'orientation',3,...
		'mat-name','xobs',...
		'buffer','1000000',...
		'position',[710,265,760,285])

add_block('built-in/Sum',[sys,'/','Sum5'])
set_param([sys,'/','Sum5'],...
		'orientation',2,...
		'inputs','++++',...
		'position',[220,139,240,221])

add_block('built-in/Sum',[sys,'/','Sum4'])
set_param([sys,'/','Sum4'],...
		'inputs','-+',...
		'position',[515,117,535,153])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain8']])
set_param([sys,'/',['Matrix',13,'Gain8']],...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:',...
		'Mask Translate','K = @1;')
set_param([sys,'/',['Matrix',13,'Gain8']],...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','F\/',...
		'position',[565,200,595,230])

add_block('built-in/MATLAB Fcn',[sys,'/',['Unv',13,'']])
set_param([sys,'/',['Unv',13,'']],...
		'MATLAB Fcn','u/(norm(u)+delta)',...
		'Mask Display','UnV',...
		'Mask Type','unit vector',...
		'Mask Dialogue','Unit Vector|delta',...
		'Mask Translate','delta=@1;',...
		'Mask Entries','0.01\/')
set_param([sys,'/',['Unv',13,'']],...
		'position',[615,200,650,230])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain7']])
set_param([sys,'/',['Matrix',13,'Gain7']],...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:',...
		'Mask Translate','K = @1;')
set_param([sys,'/',['Matrix',13,'Gain7']],...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','-rho*eye(m)\/',...
		'position',[665,200,695,230])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain6']])
set_param([sys,'/',['Matrix',13,'Gain6']],...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:',...
		'Mask Translate','K = @1;')
set_param([sys,'/',['Matrix',13,'Gain6']],...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','B\/',...
		'position',[715,138,745,162])

add_block('built-in/Sum',[sys,'/','Sum1'])
set_param([sys,'/','Sum1'],...
		'inputs','++++',...
		'position',[765,98,785,157])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain2']])
set_param([sys,'/',['Matrix',13,'Gain2']],...
		'orientation',2,...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:')
set_param([sys,'/',['Matrix',13,'Gain2']],...
		'Mask Translate','K = @1;',...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','A\/',...
		'position',[810,50,840,80])

add_block('built-in/Integrator',[sys,'/','Integrator'])
set_param([sys,'/','Integrator'],...
		'Initial','z0',...
		'position',[820,120,840,140])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain5']])
set_param([sys,'/',['Matrix',13,'Gain5']],...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:',...
		'Mask Translate','K = @1;')
set_param([sys,'/',['Matrix',13,'Gain5']],...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','-G\/',...
		'position',[600,120,630,150])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain10']])
set_param([sys,'/',['Matrix',13,'Gain10']],...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:',...
		'Mask Translate','K = @1;')
set_param([sys,'/',['Matrix',13,'Gain10']],...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','Lr\/',...
		'position',[225,245,255,275])

add_block('built-in/Mux',[sys,'/','Mux'])
set_param([sys,'/','Mux'],...
		'inputs','2',...
		'position',[150,291,180,324])

add_block('built-in/Step Fcn',[sys,'/','Step Input2'])
set_param([sys,'/','Step Input2'],...
		'Time','0',...
		'After','0.25',...
		'position',[65,270,85,290])

add_block('built-in/Step Fcn',[sys,'/','Step Input'])
set_param([sys,'/','Step Input'],...
		'Time','0',...
		'After','0',...
		'position',[65,325,85,345])

add_block('built-in/To Workspace',[sys,'/','To Workspace9'])
set_param([sys,'/','To Workspace9'],...
		'orientation',2,...
		'mat-name','u',...
		'buffer','1000000',...
		'position',[90,120,140,140])

add_block('built-in/Clock',[sys,'/','Clock'])
set_param([sys,'/','Clock'],...
		'Mask Display','',...
		'position',[40,60,60,80])

add_block('built-in/To Workspace',[sys,'/','To Workspace1'])
set_param([sys,'/','To Workspace1'],...
		'mat-name','t',...
		'buffer','1000000',...
		'position',[100,62,150,78])

add_block('built-in/State-Space',[sys,'/',['Matrix',13,'Gain9']])
set_param([sys,'/',['Matrix',13,'Gain9']],...
		'orientation',2,...
		'A','[]',...
		'B','[]',...
		'C','[]',...
		'D','K',...
		'Mask Display','K',...
		'Mask Type','Matrix Gain',...
		'Mask Dialogue','Matrix Gain.|Gain matrix:')
set_param([sys,'/',['Matrix',13,'Gain9']],...
		'Mask Translate','K = @1;',...
		'Mask Help','Multiplies input vector by entered matrix to produce output vector (y=Au).',...
		'Mask Entries','-rho*norm(S*B)*inv(S*B)\/')
set_param([sys,'/',['Matrix',13,'Gain9']],...
		'position',[280,175,310,205])
add_line(sys,[380,310;380,299;340,299;340,285])
add_line(sys,[380,299;415,299;415,240])
add_line(sys,[360,100;395,100])
add_line(sys,[375,100;375,75;480,75])
add_line(sys,[330,345;370,335])
add_line(sys,[360,345;360,370;425,370])
add_line(sys,[430,150;245,150])
add_line(sys,[845,130;890,130;890,65;845,65])
add_line(sys,[890,130;870,130;870,300;695,300])
add_line(sys,[790,300;790,345;385,335])
add_line(sys,[480,345;470,150])
add_line(sys,[415,200;415,170;245,170])
add_line(sys,[275,190;245,190])
add_line(sys,[260,260;270,260;270,210;245,210])
add_line(sys,[215,180;210,180;210,100;275,100])
add_line(sys,[210,130;145,130])
add_line(sys,[435,100;455,100;455,125;510,125])
add_line(sys,[340,195;340,190;315,190])
add_line(sys,[340,245;340,235])
add_line(sys,[655,300;500,300;510,145])
add_line(sys,[545,300;545,290])
add_line(sys,[735,300;735,290])
add_line(sys,[540,135;595,135])
add_line(sys,[550,135;560,215])
add_line(sys,[700,215;710,150])
add_line(sys,[600,215;610,215])
add_line(sys,[655,215;660,215])
add_line(sys,[230,100;230,60;575,60])
add_line(sys,[750,150;760,150])
add_line(sys,[635,135;760,135])
add_line(sys,[615,60;645,60;645,120;760,120])
add_line(sys,[805,65;740,65;740,105;760,105])
add_line(sys,[790,130;815,130])
add_line(sys,[65,70;95,70])
add_line(sys,[185,310;207,310;207,345;240,345])
add_line(sys,[207,310;207,260;220,260])
add_line(sys,[90,280;120,280;120,300;145,300])
add_line(sys,[90,335;120,335;120,315;145,315])

drawnow

% Return any arguments.
if (nargin | nargout)
	% Must use feval here to access system in memory
	if (nargin > 3)
		if (flag == 0)
			eval(['[ret,x0,str,ts,xts]=',sys,'(t,x,u,flag);'])
		else
			eval(['ret =', sys,'(t,x,u,flag);'])
		end
	else
		[ret,x0,str,ts,xts] = feval(sys);
	end
else
	drawnow % Flash up the model and execute load callback
end