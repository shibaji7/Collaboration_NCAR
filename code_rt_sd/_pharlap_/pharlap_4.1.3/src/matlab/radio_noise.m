function [out] = radio_noise(ISEASON, RLATD, ALONGD, LOCNAME, FREQA, RLMT, do_text_file, apply_correction)
%radio_noise
%
% Calculate median atmospheric noise energy at a particular location and time.
%
% ISEASON:          season index (1:winter, 2:spring, 3:summer, 4:autumn)
% RLATD:            loatitude in deg
% ALONGD:           longitude in deg
% LOCNAME:          location string
% FREQA:            model frequency in MHz
% RLMT:             local time in hours
% do_text_file:     type results to a text file 'NOISBW.OUT'
% apply_correction: logical, correction to southern hemisphere NOISEDAT
%                   software bug
%
% Author
% Mike Turley, June 2011, Defence Science & Technology Organisation
% Based on ITU-R SG3 NOISBW Recommendation ITU-R P.372
%
% C.....A FORTRAN V ROUTINE FOR PC USING USER INPUT OF SEASON, LOCATION,
% C.....FREQUENCY, BANDWIDTH AND LOCAL TIME AND 8 EXTERNAL FILES OF
% COEFFICIENTS
% C.....TO GIVE OUTPUT OF ATMOSPHERIC NOISE WITH STATISTICAL VALUES FOR DL, DU, 
% C.....SL, SM, SU, VD, AND SVD.  THIS IS A MODIFICATION (1/89) OF THE
% C.....PROGRAM NOISB.

  global ATNU ATNY CC TM RCNSE DU DL SIGM SIGU SIGL KJ JK
  global ATMO GNOS ZCNSE XADJN XNOISE ZNOISE
  global DUD FAM FAKP FAKABP
  global VDARRAY
  global DLA DUA SLA SMA SUA VD200 SVD200
  
  if nargin < 8
    apply_correction = false;
  end
  
  if nargin < 7
    do_text_file = false;
  end
  
  ATNY=0;CC=0;TM=0;RCNSE=0;DU=0;DL=0;SIGM=0;SIGU=0;SIGL=0;KJ=0;JK=0;
  DUD = zeros(5,12,5);
  FAM = zeros(14,12);
  FAKP = zeros(29,16,6);
  FAKABP = zeros(2,6);
  VDARRAY = zeros(5,12,2);
  
  SEASON = ['WINTER';'SPRING';'SUMMER';'AUTUMN'];
  SEAFIN = ['ISCOF1';'ISCOF2';'ISCOF3';'ISCOF4'];
  VDFIN = ['VDCOF1';'VDCOF2';'VDCOF3';'VDCOF4'];
  
  datapath = repmat([getenv('DIR_MODELS_REF_DAT') '/ccir_noise/'], 4, 1);
  SEAFIN = [datapath SEAFIN];
  VDFIN = [datapath VDFIN];
  
  %   TIMEBLK = ['0000-0400';'0400-0800';'0800-1200';'1200-1600';'1600-2000';'2000-2400'];
  %   TBHR = [2.0,6.0,10.0,14.0,18.0,22.0].';
  %   FREQL = [.01,.02,.05,.1,.2,.5,1.,2.,5.,10.,20.].';

  % C.....BEGIN INPUT WITH SEASON
  if nargin < 1
    ISEASON = 1; % WINTER
  end
  
  if (ISEASON == 5)
    ISE1=1;
    ISE2=4;
  else
    ISE1 = ISEASON;
    ISE2 = ISEASON;    
  end

  if nargin < 2
    RLATD = 50;  %  INPUT LOCATION LATITUDE (- IF S);
  end
  if nargin < 3
    ALONGD = 4;  %  INPUT LOCATION LONGITUDE (- IF WEST);
  end
  if nargin < 4
    LOCNAME = 'Somewhere';  % LOCATION NAME IN SINGLE QUOTES;
  end

  RLONGD = ALONGD;
  if (ALONGD < 0.0)
    RLONGD=360. + ALONGD;
  end

  % C.....INPUT FMHZ
  if nargin < 5
    FREQA = [2:2:30]; % MHz
  end
  NF = length(FREQA);

  % Local time (HH)
  if nargin < 6
    RLMT = 12.0; % Local time
  end
  
  % C.....INPUT BANDWIDTH
  BW = 200;
  IODBWDB = 1; % units: 1=dBW, 2=FA
  BWR = BW / 200.;

  % Create output file
  if do_text_file
    LUFO = fopen('NOISBW.OUT', 'wt');
  end

  % Loop over season
  for ISE = ISE1:ISE2
    ISEASON = ISE;
    KODESEA = ISEASON;
    if (RLATD < 0.0)
      if (ISEASON == 1), KODESEA = 3; end
      if (ISEASON == 2), KODESEA = 4; end
      if (ISEASON == 3), KODESEA = 1; end
      if (ISEASON == 4), KODESEA = 2; end
    end
    
    fid_3 = fopen(SEAFIN(KODESEA,:), 'r');
    temp = fread (fid_3, 'float32');
    num = numel(DUD); lo = 1; hi = lo - 1 + num;
    DUD = reshape(temp(lo:hi), [5,12,5]);
    num = numel(FAM); lo = hi+1; hi = lo-1+num;
    FAM = reshape(temp(lo:hi), [14,12]);
    num = numel(FAKP); lo = hi+1; hi = lo-1+num;
    FAKP = reshape(temp(lo:hi), [29,16,6]);
    num = numel(FAKABP); lo = hi+1; hi = lo-1+num;
    FAKABP = reshape(temp(lo:hi), [2,6]);
    fclose(fid_3);

    fid_5 = fopen(VDFIN(KODESEA,:), 'r'); 
    temp = fread (fid_5, 'float32');
    VDARRAY = reshape(temp, [5,12,2]);
    fclose(fid_5);

    % C.....USE THIS BRANCH TO DO SPECIFIC TIME
    ANOIS1(RLMT, RLATD, RLONGD);
    if do_text_file
      fprintf(LUFO,'  LAT = %6.2f,  LONG = %7.2f, %s \n', RLATD, ALONGD, LOCNAME);
      fprintf(LUFO,'  %s,  LMT = %4.1f,  BANDWIDTH = %7.0f \n', SEASON(ISEASON,:), RLMT, BW);
      if (IODBWDB == 1)
        fprintf(LUFO,'             ---MEDIAN ATMOSPHERIC NOISE OR DBW(1HZ)--\n');
      elseif (IODBWDB == 2)
        fprintf(LUFO,'             --MEDIAN ATMOSPHERIC NOISE, FA (DB>KTO)--\n');
      end
      fprintf(LUFO,'       FMHZ   NOISE   DL   DU   SL   SM   SU   VD  SVD \n');
    end
    for IFREQ = 1:NF
      FREQ = FREQA(IFREQ);
      MGENOIS(FREQ, RLATD, apply_correction)
      if (IODBWDB == 2)
        ATMO = ATMO + 204.;
      end
      VD = VDC(VD200, BWR);
      SVD = SVDC(VD200, SVD200, BWR);
      if do_text_file
        fprintf(LUFO,'    %7.3f%8.1f%5.1f%5.1f%5.1f%5.1f%5.1f%5.1f%5.1f\n', ...
          FREQ,ATMO, DLA, DUA, SLA, SMA, SUA, VD, SVD);
      end
      out.freqs(IFREQ,1) = FREQ;
      out.atmos(IFREQ,1) = ATMO;
      out.dlas(IFREQ,1) = DLA;
      out.duas(IFREQ,1) = DUA;
      out.slas(IFREQ,1) = SLA;
      out.smas(IFREQ,1) = SMA;
      out.suas(IFREQ,1) = SUA;
      out.vds(IFREQ,1) = VD;
      out.svds(IFREQ,1) = SVD;
    end

  end % Season

  if do_text_file
    fclose(LUFO);
  end

  % Output
  out.season_id = ISEASON;
  out.season_name = SEASON(ISEASON,:);
  out.lat = RLATD;
  out.lon = RLONGD;
  out.location = LOCNAME;
  out.local_time = RLMT;
  out.bandwidth_hz = BW;
  noise_types = {'dBW', 'FA'};
  out.noise_type = char(noise_types{IODBWDB}); % units: 1=dBW, 2=FA
  
return
      
function [] = ANOIS1(RLMT, RLATD, RLONGD)
  % CR....A ROUTINE THAT USES RLMT TO DETERMINE THE TIMEBLOCK (KJ)
  % CR....AND ADJACENT TIME BLOCK (JK) (THIS IS THE PRIOR TIMEBLOCK
  % CR....FOR THE FIRST 2 HOURS OF KJ, THE SAME, IE JK=KJ, FOR THE 3RD
  % CR....HOUR OF KJ AND THE NEXT TIME BLOCK FOR THE LAST HOUR OF KJ)
  % CR....AND THEN CALLS NOISY TO FIGURE THE ATMOSPHERIC NOISE (ATNU
  % CR....OR ATNY) FOR EACH OF THESE TIME BLOCKS.
  % C.....
  % C.....THIS ROUTINE DETERMINES THE 1 MHZ ATMOSPHERIC NOISE
  % C.....
  % C.....FOURIER SERIES IN LATITUDE AND LONGITUDE FOR TWO DISCRETE
  % C.....LOCAL TIME BLOCKS
  % C.....
  
  global ATNU ATNY CC TM RCNSE DU DL SIGM SIGU SIGL KJ JK
  
  % C.....LMT AT RCVR SITE
  CC = RLMT;
  KJ = 6;
  if (CC - 22.) < 0
    KJ = fix(CC / 4. + 1.);
  end

  TM = 4 * KJ - 2;

  if (CC - TM) < 0
    JK = KJ - 1;
  elseif (CC - TM) == 0
    JK = KJ;
  else
    JK = KJ + 1;
  end
  
  if (JK <= 0)
    JK = 6;
  else
    if (JK - 6) > 0
      JK = 1;
    end
  end
  % C.....EAST LONGITUDE (IN DEGREES)
  CEG = RLONGD;

  XLA = RLATD;
  % C.....LATITUDE (IN DEGREES) "+" IS NORTH
  ATNU = NOISY(KJ, XLA, CEG);
  ATNY = NOISY(JK, XLA, CEG);
return 

function [ANOS] = NOISY (KJ, XLA, CEG) 
  % CR....A ROUTINE TO USE THE TIMEBLOCK (KJ), THE LAT (XLA), THE LONG
  % CR....(CEG), AND THE COEFFICIENTS (FAKP AND FAKAB) TO
  % CR....DETERMINE THE ATMOSPHERIC NOISE (ANOS).
  % CR....THIS ROUTINE USES MAPS TO GET THE 1 MHZ FAM VALUE.
  % C.....NOISY IS A GENERAL PURPOSE ROUTINE USED TO EVALUATE A FOURIER
  % C.....SERIES IN TWO VARIABLES.
  % C.....KJ --- NUMBER OF FOURIER COEFFICIENT ARRAY TO BE USED
  % C.....XLA --- GEOGRAPHIC LATITUDE, DEGREES,
  % C.....CEG --- GEOGRAPHIC EAST LONGITUDE, DEGREES
  % C.....ANOS --- NOISE VALUE, MEDIAN POWER DB ABOVE KTB
  % C.....FAKABP --- NORMALIZING FACTORS FOR FOURIER SERIES
  % C.....KJ = 1 TO 6 IS ATMOSPHERIC NOISE
  % C.....
  % C.....* NOTE - XLA, CEG, ANOS, FAKABP ARE NOT ALWAYS AS PREVIOUSLY
  % C.....         DEFINED
  % C.....FOURIER VARIABLES AND ATMOSPHERIC RADIO NOISE
  % C.....
  
  global DUD FAM FAKP FAKABP
  
  SX = zeros(15,1);
  SY = zeros(29,1);
  ZZ = zeros(29,1);
  if (KJ - 6) > 0
    KJ = 6;
  end
  % C.....LIMITS OF FOURIER SERIES
  LM = 29;
  LN = 15;
  % C.....HALF ANGLE (IN RADIANS)
  Q = .0087266466 * CEG;

  % C.....LONGITUDE SINES
  SX = sin(Q * [1:LN]);

  % C.....LONGITUDE SERIES
  ZZ = FAKP(:,1:16,KJ) * [SX, 1].';
  % C.....ANGLE PLUS 90 DEGREES (IN RADIANS)
  Q = .01745329252 * (XLA + 90.);
  % C.....LATITUDE SERIES
  SY = sin(Q * [1:29]);
  R = SY(1:LM) * ZZ(1:LM);
  
  % C.....FINAL FOURIER SERIES EVALUATION (NOTE LINEAR NORMALIZATION)
  ANOS = R + FAKABP(1,KJ) + FAKABP(2,KJ)* Q;
return
      
function [] = MGENOIS(FREQ, RLAT, apply_correction)
  % C.....
  % CR....THIS ROUTINE IS A MODIFIED VERSION OF GENOIS WHERE JUST THE ATMOSPHERIC
  % CR....NOISE STATISTICS ARE CALCULATED AND VD AND SVD ARE ADDED.
  % C.....
  % %       COMMON /ANOIS/ ATNU,ATNY,CC,TM,RCNSE,DU,DL,SIGM,SIGU,SIGL,KJ,JK
  %       COMMON /TON/ ATMO, GNOS, ZCNSE, XADJN, XNOISE, ZNOISE
  %       COMMON /NSTAT/ DLA,DUA,SLA,SMA,SUA,VD200,SVD200
  global ATNU ATNY CC TM RCNSE DU DL SIGM SIGU SIGL KJ JK
  global ATMO GNOS ZCNSE XADJN XNOISE ZNOISE
  global DLA DUA SLA SMA SUA VD200 SVD200
  
  DUME = min([FREQ, 55.]);
  [ATNZ,DU,DL,SIGM,SIGU,SIGL,VD1,SVD1] = MGENFAM(RLAT,KJ,DUME,ATNU, apply_correction);
  [ATNX,DX,DQ,SIGZ,SIGX,SIGSQ,VD2,SVD2] = MGENFAM(RLAT,JK,DUME,ATNY, apply_correction);
  
  % C.....BEGIN OF INTERPOLATION ON LOCAL TIME
  SLOP = abs(CC-TM)/4.;
  ATNOS = ATNZ + (ATNX - ATNZ) * SLOP;
  ATMO = ATNOS - 204.;
  DUA = DU + (DX-DU)*SLOP;
  DLA = DL + (DQ-DL)*SLOP;
  SMA = SIGM + (SIGZ-SIGM)*SLOP;
  SUA = SIGU + (SIGX-SIGU)*SLOP;
  SLA = SIGL + (SIGSQ-SIGL)* SLOP;
  VD200 = VD1 + (VD2-VD1)*SLOP;
  SVD200 = SVD1 + (SVD2-SVD1)*SLOP;
return
      
function [FA,DU,DL,DMS,DUS,DLS,VD,SVD] = MGENFAM(Y2,IBLK,FREQ,Z, apply_correction)
  % c**********************************************************************
  % c          Re-written 3.June.93 by Greg Hand because previous version was
  % c          really incorrect. It made an attempt to limit Sigma Fam (DMS)
  % c          to a 10 MHz frequency, but the indicies I and J became
  % c          confused, and the result was not correct. This current
  % c          version should limit DMS to 10 MHz and the others to 20 MHz
  % c          because the curves end at 20 MHz. Unfortunately, this error
  % c          has probably existed since time began, and it may take a
  % c          while for this corrected version to propagate into all version
  % c          that exist. The magnitude of the error that would have been
  % c          caused is not known, but it is believed to be small.
  % C
  % C Amendment
  % C          Turley June 2011: apply_correction to southern hemisphere calculations 
  % C
  % c**********************************************************************
  % C.....
  % CR....THIS IS A MODIFIED GENFAM ROUTINE; AN ADDITIONAL SECTION DEFINES VD, SVD.
  % C.....GENFAM CALCULATES THE FREQUENCY DEPENDENCE OF THE ATMOSPHERIC
  % C.....NOISE AND GETS DECILES AND PREDICTION ERRORS FROM TABLES
  % C.....
  %   COMMON /TWO/ DUD(5,12,5),FAM(14,12),FAKP(29,16,6),FAKABP(2,6)
  %   COMMON /THREE/ VDARRAY(5,12,2)
  
  global DUD FAM FAKP FAKABP
  global VDARRAY
  
  V = zeros(5,1);
  IBK = IBLK;
  % C.....CHECK IF LATITUDE IS NORTH OR SOUTH
  %MDET corrects what appeared to index the wrong seasonal dependent freq conversion coefficients
  if ~apply_correction
    if (Y2 < 0.)
      IBK = IBK + 6;
    end
  end
  %MDET
  U1 = - .75;
  X = log10(FREQ);
  U = (8. * 2.^X - 11.) / 4.;
  for KOP = 1:2
    PZ = U1 * FAM (1, IBK) + FAM (2, IBK);
    PX = U1 * FAM (8, IBK) + FAM (9, IBK);
    for I = 3:7
      PZ = U1 * PZ + FAM (I, IBK);
      PX = U1 * PX + FAM (I + 7, IBK);
    end
    if (KOP == 1)
      CZ = Z * PZ + PX;
      CZ = Z + Z - CZ;
      U1 = U;
    end
  end
  FA = CZ * PZ + PX;
  % c          Limit frequency to 20 MHz for DUA, DLA, DUS, DLS
  % c            because curves in REP 322 only go to 20 MHz
  if(FREQ > 20.)
    X = log10(20.);
  end
  x_mat = repmat(X.^(4:-1:0).', [1,5]);
  % c          Limit frequency to 10 MHz for DMS (Sigma Fam)
  % c            because curves in REP 322 only go to 10 MHz
  if FREQ > 10.
    x_mat(:,5) = 1;
  end
  V = sum(x_mat .* squeeze(DUD(1:5, IBK, 1:5)));
  DU = V(1);
  DL = V(2);
  DUS = V(3);
  DLS = V(4);
  DMS = V(5);
  X = log10(FREQ);
  if(FREQ > 20.)
    X = log10(20.);
  end
  x_vec = X.^(4:-1:0);
  V = x_vec * squeeze(VDARRAY (1:5, IBK, 1:2));
  VD = V(1);
  SVD = V(2);
return

function [VDC1] = VDC(VD200, BWR)
  % C     OBTAINS THE NOISE PARAMETER -VD-, FOR THE SPECIFIED BANDWIDTH
  % C     FROM THE CCIR REPORT 322 (OR OTHER) 200HZ BANDWIDTH -VD- (VD200),
  % C     -BWR- IS THE BANDWIDTH RATIO (REQUIRED BANDWIDTH/200 HZ BANDWIDTH).
  VDC1 = 1.049;
  if (VD200 <= 1.049)
    return
  end
  VDO = VD200 + (0.4679 + 0.2111 * VD200) * log10(BWR);
  if (VDO <= 1.049)
    return
  end
  VDC1 = VDO;
return
    
function [SVDC1] = SVDC(VD200, SVD200, BWR)
  V01 = VDC((VD200 + SVD200), BWR);
  V02 = VDC((VD200 - SVD200), BWR);
  SVDC1 = (V01 - V02) / 2.;
return

