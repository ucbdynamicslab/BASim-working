{
double __53844,__53845,__53846,__53847,__53848,__53849,__53852,__53853,__53854,__53855;

__53844 = -yj;
__53845 = yi + __53844;
__53846 = Power(__53845,2);
__53847 = -zj;
__53848 = zi + __53847;
__53849 = Power(__53848,2);
__53852 = -xj;
__53853 = xi + __53852;
__53854 = Power(__53853,2);
__53855 = __53854 + __53846 + __53849;
mathematicaJacobian(0,0) = (ks*(len*(__53846 + __53849) - Power(__53855,1.5)))/(len*Power(__53855,1.5));
}

{
double __53862,__53863,__53864,__53865;

__53862 = -xj;
__53863 = xi + __53862;
__53864 = -yj;
__53865 = yi + __53864;
mathematicaJacobian(0,1) = -((ks*__53863*__53865)/Power(Power(__53863,2) + Power(__53865,2) + Power(zi - zj,2),1.5));
}

{
double __53875,__53876,__53881,__53882;

__53875 = -xj;
__53876 = xi + __53875;
__53881 = -zj;
__53882 = zi + __53881;
mathematicaJacobian(0,2) = -((ks*__53876*__53882)/Power(Power(__53876,2) + Power(yi - yj,2) + Power(__53882,2),1.5));
}

{
double __53889,__53890,__53891,__53892,__53893,__53894,__53897,__53898,__53899,__53900;

__53889 = -yj;
__53890 = yi + __53889;
__53891 = Power(__53890,2);
__53892 = -zj;
__53893 = zi + __53892;
__53894 = Power(__53893,2);
__53897 = -xj;
__53898 = xi + __53897;
__53899 = Power(__53898,2);
__53900 = __53899 + __53891 + __53894;
mathematicaJacobian(0,3) = -((ks*(len*(__53891 + __53894) - Power(__53900,1.5)))/(len*Power(__53900,1.5)));
}

{
double __53907,__53908,__53909,__53910;

__53907 = -xj;
__53908 = xi + __53907;
__53909 = -yj;
__53910 = yi + __53909;
mathematicaJacobian(0,4) = (ks*__53908*__53910)/Power(Power(__53908,2) + Power(__53910,2) + Power(zi - zj,2),1.5);
}

{
double __53920,__53921,__53926,__53927;

__53920 = -xj;
__53921 = xi + __53920;
__53926 = -zj;
__53927 = zi + __53926;
mathematicaJacobian(0,5) = (ks*__53921*__53927)/Power(Power(__53921,2) + Power(yi - yj,2) + Power(__53927,2),1.5);
}

{
double __53933,__53934,__53935,__53936;

__53933 = -xj;
__53934 = xi + __53933;
__53935 = -yj;
__53936 = yi + __53935;
mathematicaJacobian(1,0) = -((ks*__53934*__53936)/Power(Power(__53934,2) + Power(__53936,2) + Power(zi - zj,2),1.5));
}

{
double __53947,__53948,__53949,__53950,__53951,__53952,__53955,__53956,__53957,__53958;

__53947 = -xj;
__53948 = xi + __53947;
__53949 = Power(__53948,2);
__53950 = -zj;
__53951 = zi + __53950;
__53952 = Power(__53951,2);
__53955 = -yj;
__53956 = yi + __53955;
__53957 = Power(__53956,2);
__53958 = __53949 + __53957 + __53952;
mathematicaJacobian(1,1) = (ks*(len*(__53949 + __53952) - Power(__53958,1.5)))/(len*Power(__53958,1.5));
}

{
double __53965,__53966,__53971,__53972;

__53965 = -yj;
__53966 = yi + __53965;
__53971 = -zj;
__53972 = zi + __53971;
mathematicaJacobian(1,2) = -((ks*__53966*__53972)/Power(Power(xi - xj,2) + Power(__53966,2) + Power(__53972,2),1.5));
}

{
double __53978,__53979,__53980,__53981;

__53978 = -xj;
__53979 = xi + __53978;
__53980 = -yj;
__53981 = yi + __53980;
mathematicaJacobian(1,3) = (ks*__53979*__53981)/Power(Power(__53979,2) + Power(__53981,2) + Power(zi - zj,2),1.5);
}

{
double __53992,__53993,__53994,__53995,__53996,__53997,__54000,__54001,__54002,__54003;

__53992 = -xj;
__53993 = xi + __53992;
__53994 = Power(__53993,2);
__53995 = -zj;
__53996 = zi + __53995;
__53997 = Power(__53996,2);
__54000 = -yj;
__54001 = yi + __54000;
__54002 = Power(__54001,2);
__54003 = __53994 + __54002 + __53997;
mathematicaJacobian(1,4) = -((ks*(len*(__53994 + __53997) - Power(__54003,1.5)))/(len*Power(__54003,1.5)));
}

{
double __54010,__54011,__54016,__54017;

__54010 = -yj;
__54011 = yi + __54010;
__54016 = -zj;
__54017 = zi + __54016;
mathematicaJacobian(1,5) = (ks*__54011*__54017)/Power(Power(xi - xj,2) + Power(__54011,2) + Power(__54017,2),1.5);
}

{
double __54023,__54024,__54029,__54030;

__54023 = -xj;
__54024 = xi + __54023;
__54029 = -zj;
__54030 = zi + __54029;
mathematicaJacobian(2,0) = -((ks*__54024*__54030)/Power(Power(__54024,2) + Power(yi - yj,2) + Power(__54030,2),1.5));
}

{
double __54036,__54037,__54042,__54043;

__54036 = -yj;
__54037 = yi + __54036;
__54042 = -zj;
__54043 = zi + __54042;
mathematicaJacobian(2,1) = -((ks*__54037*__54043)/Power(Power(xi - xj,2) + Power(__54037,2) + Power(__54043,2),1.5));
}

{
double __54050,__54051,__54052,__54053,__54054,__54055,__54058,__54059,__54060,__54061;

__54050 = -xj;
__54051 = xi + __54050;
__54052 = Power(__54051,2);
__54053 = -yj;
__54054 = yi + __54053;
__54055 = Power(__54054,2);
__54058 = -zj;
__54059 = zi + __54058;
__54060 = Power(__54059,2);
__54061 = __54052 + __54055 + __54060;
mathematicaJacobian(2,2) = (ks*(len*(__54052 + __54055) - Power(__54061,1.5)))/(len*Power(__54061,1.5));
}

{
double __54068,__54069,__54074,__54075;

__54068 = -xj;
__54069 = xi + __54068;
__54074 = -zj;
__54075 = zi + __54074;
mathematicaJacobian(2,3) = (ks*__54069*__54075)/Power(Power(__54069,2) + Power(yi - yj,2) + Power(__54075,2),1.5);
}

{
double __54081,__54082,__54087,__54088;

__54081 = -yj;
__54082 = yi + __54081;
__54087 = -zj;
__54088 = zi + __54087;
mathematicaJacobian(2,4) = (ks*__54082*__54088)/Power(Power(xi - xj,2) + Power(__54082,2) + Power(__54088,2),1.5);
}

{
double __54095,__54096,__54097,__54098,__54099,__54100,__54103,__54104,__54105,__54106;

__54095 = -xj;
__54096 = xi + __54095;
__54097 = Power(__54096,2);
__54098 = -yj;
__54099 = yi + __54098;
__54100 = Power(__54099,2);
__54103 = -zj;
__54104 = zi + __54103;
__54105 = Power(__54104,2);
__54106 = __54097 + __54100 + __54105;
mathematicaJacobian(2,5) = -((ks*(len*(__54097 + __54100) - Power(__54106,1.5)))/(len*Power(__54106,1.5)));
}

{
double __54114,__54115,__54116,__54117,__54118,__54119,__54122,__54123,__54124,__54125;

__54114 = -yj;
__54115 = yi + __54114;
__54116 = Power(__54115,2);
__54117 = -zj;
__54118 = zi + __54117;
__54119 = Power(__54118,2);
__54122 = -xj;
__54123 = xi + __54122;
__54124 = Power(__54123,2);
__54125 = __54124 + __54116 + __54119;
mathematicaJacobian(3,0) = -((ks*(len*(__54116 + __54119) - Power(__54125,1.5)))/(len*Power(__54125,1.5)));
}

{
double __54132,__54133,__54134,__54135;

__54132 = -xj;
__54133 = xi + __54132;
__54134 = -yj;
__54135 = yi + __54134;
mathematicaJacobian(3,1) = (ks*__54133*__54135)/Power(Power(__54133,2) + Power(__54135,2) + Power(zi - zj,2),1.5);
}

{
double __54145,__54146,__54151,__54152;

__54145 = -xj;
__54146 = xi + __54145;
__54151 = -zj;
__54152 = zi + __54151;
mathematicaJacobian(3,2) = (ks*__54146*__54152)/Power(Power(__54146,2) + Power(yi - yj,2) + Power(__54152,2),1.5);
}

{
double __54159,__54160,__54161,__54162,__54163,__54164,__54167,__54168,__54169,__54170;

__54159 = -yj;
__54160 = yi + __54159;
__54161 = Power(__54160,2);
__54162 = -zj;
__54163 = zi + __54162;
__54164 = Power(__54163,2);
__54167 = -xj;
__54168 = xi + __54167;
__54169 = Power(__54168,2);
__54170 = __54169 + __54161 + __54164;
mathematicaJacobian(3,3) = (ks*(len*(__54161 + __54164) - Power(__54170,1.5)))/(len*Power(__54170,1.5));
}

{
double __54177,__54178,__54179,__54180;

__54177 = -xj;
__54178 = xi + __54177;
__54179 = -yj;
__54180 = yi + __54179;
mathematicaJacobian(3,4) = -((ks*__54178*__54180)/Power(Power(__54178,2) + Power(__54180,2) + Power(zi - zj,2),1.5));
}

{
double __54190,__54191,__54196,__54197;

__54190 = -xj;
__54191 = xi + __54190;
__54196 = -zj;
__54197 = zi + __54196;
mathematicaJacobian(3,5) = -((ks*__54191*__54197)/Power(Power(__54191,2) + Power(yi - yj,2) + Power(__54197,2),1.5));
}

{
double __54203,__54204,__54205,__54206;

__54203 = -xj;
__54204 = xi + __54203;
__54205 = -yj;
__54206 = yi + __54205;
mathematicaJacobian(4,0) = (ks*__54204*__54206)/Power(Power(__54204,2) + Power(__54206,2) + Power(zi - zj,2),1.5);
}

{
double __54217,__54218,__54219,__54220,__54221,__54222,__54225,__54226,__54227,__54228;

__54217 = -xj;
__54218 = xi + __54217;
__54219 = Power(__54218,2);
__54220 = -zj;
__54221 = zi + __54220;
__54222 = Power(__54221,2);
__54225 = -yj;
__54226 = yi + __54225;
__54227 = Power(__54226,2);
__54228 = __54219 + __54227 + __54222;
mathematicaJacobian(4,1) = -((ks*(len*(__54219 + __54222) - Power(__54228,1.5)))/(len*Power(__54228,1.5)));
}

{
double __54235,__54236,__54241,__54242;

__54235 = -yj;
__54236 = yi + __54235;
__54241 = -zj;
__54242 = zi + __54241;
mathematicaJacobian(4,2) = (ks*__54236*__54242)/Power(Power(xi - xj,2) + Power(__54236,2) + Power(__54242,2),1.5);
}

{
double __54248,__54249,__54250,__54251;

__54248 = -xj;
__54249 = xi + __54248;
__54250 = -yj;
__54251 = yi + __54250;
mathematicaJacobian(4,3) = -((ks*__54249*__54251)/Power(Power(__54249,2) + Power(__54251,2) + Power(zi - zj,2),1.5));
}

{
double __54262,__54263,__54264,__54265,__54266,__54267,__54270,__54271,__54272,__54273;

__54262 = -xj;
__54263 = xi + __54262;
__54264 = Power(__54263,2);
__54265 = -zj;
__54266 = zi + __54265;
__54267 = Power(__54266,2);
__54270 = -yj;
__54271 = yi + __54270;
__54272 = Power(__54271,2);
__54273 = __54264 + __54272 + __54267;
mathematicaJacobian(4,4) = (ks*(len*(__54264 + __54267) - Power(__54273,1.5)))/(len*Power(__54273,1.5));
}

{
double __54280,__54281,__54286,__54287;

__54280 = -yj;
__54281 = yi + __54280;
__54286 = -zj;
__54287 = zi + __54286;
mathematicaJacobian(4,5) = -((ks*__54281*__54287)/Power(Power(xi - xj,2) + Power(__54281,2) + Power(__54287,2),1.5));
}

{
double __54293,__54294,__54299,__54300;

__54293 = -xj;
__54294 = xi + __54293;
__54299 = -zj;
__54300 = zi + __54299;
mathematicaJacobian(5,0) = (ks*__54294*__54300)/Power(Power(__54294,2) + Power(yi - yj,2) + Power(__54300,2),1.5);
}

{
double __54306,__54307,__54312,__54313;

__54306 = -yj;
__54307 = yi + __54306;
__54312 = -zj;
__54313 = zi + __54312;
mathematicaJacobian(5,1) = (ks*__54307*__54313)/Power(Power(xi - xj,2) + Power(__54307,2) + Power(__54313,2),1.5);
}

{
double __54320,__54321,__54322,__54323,__54324,__54325,__54328,__54329,__54330,__54331;

__54320 = -xj;
__54321 = xi + __54320;
__54322 = Power(__54321,2);
__54323 = -yj;
__54324 = yi + __54323;
__54325 = Power(__54324,2);
__54328 = -zj;
__54329 = zi + __54328;
__54330 = Power(__54329,2);
__54331 = __54322 + __54325 + __54330;
mathematicaJacobian(5,2) = -((ks*(len*(__54322 + __54325) - Power(__54331,1.5)))/(len*Power(__54331,1.5)));
}

{
double __54338,__54339,__54344,__54345;

__54338 = -xj;
__54339 = xi + __54338;
__54344 = -zj;
__54345 = zi + __54344;
mathematicaJacobian(5,3) = -((ks*__54339*__54345)/Power(Power(__54339,2) + Power(yi - yj,2) + Power(__54345,2),1.5));
}

{
double __54351,__54352,__54357,__54358;

__54351 = -yj;
__54352 = yi + __54351;
__54357 = -zj;
__54358 = zi + __54357;
mathematicaJacobian(5,4) = -((ks*__54352*__54358)/Power(Power(xi - xj,2) + Power(__54352,2) + Power(__54358,2),1.5));
}

{
double __54365,__54366,__54367,__54368,__54369,__54370,__54373,__54374,__54375,__54376;

__54365 = -xj;
__54366 = xi + __54365;
__54367 = Power(__54366,2);
__54368 = -yj;
__54369 = yi + __54368;
__54370 = Power(__54369,2);
__54373 = -zj;
__54374 = zi + __54373;
__54375 = Power(__54374,2);
__54376 = __54367 + __54370 + __54375;
mathematicaJacobian(5,5) = (ks*(len*(__54367 + __54370) - Power(__54376,1.5)))/(len*Power(__54376,1.5));
}

