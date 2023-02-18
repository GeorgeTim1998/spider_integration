import numpy as np
from math import pi, sqrt
import matplotlib.pyplot as pyplot

ACASE48 = 'KIAM'
SPIDER = 'SPIDER'
IDUM = 3
MESHR = 128
MESHZ = 257

folder = 'Files'
filename = 'eqdsk_equilx'
# with open("%s/%s" % (folder, filename), 'w') as file:
#   file.write("%8s%8s%36d%4d%4d" % (ACASE48, SPIDER, IDUM, MESHR, MESHZ))

fpol = np.array([0.323515213E+01, 0.322502855E+01, 0.321626699E+01, 0.320842266E+01, 0.320102936E+01,
0.319397791E+01, 0.318720758E+01, 0.318067832E+01, 0.317435701E+01, 0.316822879E+01,
0.316227767E+01, 0.315649170E+01, 0.315086024E+01, 0.314537422E+01, 0.314002260E+01,
0.313480080E+01, 0.312970310E+01, 0.312472346E+01, 0.311985839E+01, 0.311510409E+01,
0.311045572E+01, 0.310590983E+01, 0.310146411E+01, 0.309711598E+01, 0.309286230E+01,
0.308870002E+01, 0.308462672E+01, 0.308064063E+01, 0.307674010E+01, 0.307292314E+01,
0.306918761E+01, 0.306553180E+01, 0.306195442E+01, 0.305845419E+01, 0.305502990E+01,
0.305168036E+01, 0.304840430E+01, 0.304520027E+01, 0.304206678E+01, 0.303900265E+01,
0.303600688E+01, 0.303307851E+01, 0.303021667E+01, 0.302742056E+01, 0.302468941E+01,
0.302202235E+01, 0.301941844E+01, 0.301687668E+01, 0.301439617E+01, 0.301197612E+01,
0.300961580E+01, 0.300731457E+01, 0.300507179E+01, 0.300288689E+01, 0.300075929E+01,
0.299868836E+01, 0.299667333E+01, 0.299471339E+01, 0.299280779E+01, 0.299095595E+01,
0.298915738E+01, 0.298741158E+01, 0.298571805E+01, 0.298407631E+01, 0.298248588E+01,
0.298094618E+01, 0.297945646E+01, 0.297801605E+01, 0.297662441E+01, 0.297528108E+01,
0.297398561E+01, 0.297273757E+01, 0.297153648E+01, 0.297038187E+01, 0.296927323E+01,
0.296820999E+01, 0.296719155E+01, 0.296621734E+01, 0.296528685E+01, 0.296439963E+01,
0.296355522E+01, 0.296275320E+01, 0.296199311E+01, 0.296127447E+01, 0.296059671E+01,
0.295995916E+01, 0.295936118E+01, 0.295880225E+01, 0.295828189E+01, 0.295779962E+01,
0.295735496E+01, 0.295694730E+01, 0.295657594E+01, 0.295624025E+01, 0.295593974E+01,
0.295567387E+01, 0.295544199E+01, 0.295524330E+01, 0.295507697E+01, 0.295494230E+01,
0.295483864E+01, 0.295476519E+01, 0.295472088E+01, 0.295470479E+01, 0.295471629E+01,
0.295475447E+01, 0.295481803E+01, 0.295490593E+01, 0.295501735E+01, 0.295515096E+01,
0.295530547E+01, 0.295547998E+01, 0.295567293E+01, 0.295588281E+01, 0.295610856E+01,
0.295634866E+01, 0.295660178E+01, 0.295686663E+01, 0.295714176E+01, 0.295742605E+01,
0.295771844E+01, 0.295801820E+01, 0.295832518E+01, 0.295863961E+01, 0.295896263E+01,
0.295929569E+01, 0.295964087E+01, 0.296000000E+01])

pres = np.array([0.798439125E+05, 0.788314721E+05, 0.781283908E+05, 0.775965912E+05, 0.770949890E+05,
0.766032601E+05, 0.761164706E+05, 0.756325697E+05, 0.751498775E+05, 0.746677913E+05,
0.741856883E+05, 0.737031763E+05, 0.732199141E+05, 0.727356587E+05, 0.722501362E+05,
0.717632391E+05, 0.712748489E+05, 0.707848471E+05, 0.702931688E+05, 0.697997457E+05,
0.693044993E+05, 0.688073748E+05, 0.683083353E+05, 0.678073406E+05, 0.673043454E+05,
0.667993062E+05, 0.662921877E+05, 0.657829629E+05, 0.652716063E+05, 0.647580876E+05,
0.642423732E+05, 0.637244352E+05, 0.632042499E+05, 0.626817940E+05, 0.621570444E+05,
0.616299779E+05, 0.611005694E+05, 0.605687886E+05, 0.600346031E+05, 0.594979842E+05,
0.589589054E+05, 0.584173408E+05, 0.578732651E+05, 0.573266539E+05, 0.567774826E+05,
0.562257234E+05, 0.556713446E+05, 0.551143118E+05, 0.545545902E+05, 0.539921475E+05,
0.534269527E+05, 0.528589758E+05, 0.522881874E+05, 0.517145585E+05, 0.511380601E+05,
0.505586592E+05, 0.499763160E+05, 0.493909854E+05, 0.488026232E+05, 0.482111923E+05,
0.476166594E+05, 0.470189906E+05, 0.464181509E+05, 0.458141053E+05, 0.452068186E+05,
0.445962475E+05, 0.439823368E+05, 0.433650324E+05, 0.427442889E+05, 0.421200657E+05,
0.414923230E+05, 0.408610204E+05, 0.402261153E+05, 0.395875625E+05, 0.389453136E+05,
0.382993158E+05, 0.376495106E+05, 0.369958410E+05, 0.363382557E+05, 0.356767079E+05,
0.350111527E+05, 0.343415444E+05, 0.336678382E+05, 0.329899860E+05, 0.323079299E+05,
0.316216069E+05, 0.309309592E+05, 0.302359408E+05, 0.295365139E+05, 0.288326419E+05,
0.281242888E+05, 0.274114173E+05, 0.266939918E+05, 0.259719889E+05, 0.252453964E+05,
0.245142054E+05, 0.237784118E+05, 0.230380287E+05, 0.222930924E+05, 0.215436508E+05,
0.207897529E+05, 0.200314833E+05, 0.192689891E+05, 0.185024206E+05, 0.177318963E+05,
0.169576056E+05, 0.161798618E+05, 0.153989406E+05, 0.146150893E+05, 0.138287283E+05,
0.130403106E+05, 0.122501780E+05, 0.114589613E+05, 0.106673189E+05, 0.987576351E+04,
0.908506424E+04, 0.829595271E+04, 0.750921086E+04, 0.672582581E+04, 0.594669670E+04,
0.517294738E+04, 0.440562935E+04, 0.364584623E+04, 0.289473510E+04, 0.215341297E+04,
0.142304947E+04, 0.704892223E+03, 0.000000000E+00])

ffprim = np.array([0.999417939E+01, 0.867587335E+01, 0.743508148E+01, 0.695207210E+01, 0.657569454E+01,
0.629749948E+01, 0.603575254E+01, 0.583357886E+01, 0.563452217E+01, 0.546150704E+01,
0.529334016E+01, 0.514365366E+01, 0.499626043E+01, 0.486438730E+01, 0.473811418E+01,
0.461570973E+01, 0.450077253E+01, 0.439026714E+01, 0.428251854E+01, 0.417977053E+01,
0.408147937E+01, 0.398571004E+01, 0.389229592E+01, 0.380187895E+01, 0.371470496E+01,
0.363036779E+01, 0.354810259E+01, 0.346752112E+01, 0.338872287E+01, 0.331221237E+01,
0.323773116E+01, 0.316469114E+01, 0.309300816E+01, 0.302264966E+01, 0.295356410E+01,
0.288572434E+01, 0.281925283E+01, 0.275437377E+01, 0.269086394E+01, 0.262842895E+01,
0.256701983E+01, 0.250657815E+01, 0.244701795E+01, 0.238828344E+01, 0.233038564E+01,
0.227343443E+01, 0.221747784E+01, 0.216249287E+01, 0.210837582E+01, 0.205502311E+01,
0.200237930E+01, 0.195040213E+01, 0.189905573E+01, 0.184830359E+01, 0.179815156E+01,
0.174870963E+01, 0.170005559E+01, 0.165219784E+01, 0.160500625E+01, 0.155832478E+01,
0.151214135E+01, 0.146646540E+01, 0.142129133E+01, 0.137659656E+01, 0.133238359E+01,
0.128882043E+01, 0.124599020E+01, 0.120372527E+01, 0.116192131E+01, 0.112055432E+01,
0.107961533E+01, 0.103910788E+01, 0.999052991E+00, 0.959468650E+00, 0.920373405E+00,
0.881818016E+00, 0.843823765E+00, 0.806325328E+00, 0.769264603E+00, 0.732606981E+00,
0.696357684E+00, 0.660517176E+00, 0.625070175E+00, 0.590093788E+00, 0.555689447E+00,
0.521878675E+00, 0.488592171E+00, 0.455719963E+00, 0.423260670E+00, 0.391226511E+00,
0.359638374E+00, 0.328646675E+00, 0.298274038E+00, 0.268349746E+00, 0.238847659E+00,
0.209839462E+00, 0.181451176E+00, 0.153801161E+00, 0.126834781E+00, 0.100389498E+00,
0.744817446E-01, 0.494291817E-01, 0.253368376E-01, 0.180020262E-02,-0.212091973E-01,
-0.431993615E-01, -0.640161112E-01, -0.842876889E-01, -0.103681598E+00, -0.121748344E+00,
-0.138964642E+00, -0.155487616E+00, -0.170113015E+00, -0.184089489E+00, -0.196908518E+00,
-0.208350897E+00, -0.218901438E+00, -0.228207261E+00, -0.236299558E+00, -0.243714174E+00,
-0.250018929E+00, -0.256210321E+00, -0.262269054E+00, -0.268844714E+00, -0.276835486E+00,
-0.285875810E+00, -0.297368729E+00, -0.309533860E+00])

pprime = np.array([0.335012952E+06, 0.243038167E+06, 0.158384855E+06, 0.145245202E+06, 0.141143698E+06,
0.139608090E+06, 0.138323506E+06, 0.137958905E+06, 0.137633363E+06, 0.137612867E+06,
0.137643056E+06, 0.137846315E+06, 0.138071389E+06, 0.138413402E+06, 0.138794793E+06,
0.139198297E+06, 0.139647210E+06, 0.140118480E+06, 0.140604390E+06, 0.141114657E+06,
0.141645411E+06, 0.142186965E+06, 0.142738784E+06, 0.143303236E+06, 0.143881002E+06,
0.144470289E+06, 0.145068112E+06, 0.145672922E+06, 0.146285204E+06, 0.146907424E+06,
0.147538788E+06, 0.148177022E+06, 0.148821883E+06, 0.149473384E+06, 0.150131487E+06,
0.150796217E+06, 0.151468649E+06, 0.152150683E+06, 0.152841584E+06, 0.153540067E+06,
0.154246029E+06, 0.154959339E+06, 0.155679726E+06, 0.156406983E+06, 0.157141464E+06,
0.157884486E+06, 0.158637097E+06, 0.159399845E+06, 0.160172202E+06, 0.160953504E+06,
0.161743469E+06, 0.162541926E+06, 0.163348741E+06, 0.164163671E+06, 0.164987074E+06,
0.165820882E+06, 0.166666964E+06, 0.167526508E+06, 0.168397855E+06, 0.169278629E+06,
0.170168890E+06, 0.171069065E+06, 0.171979329E+06, 0.172899457E+06, 0.173829813E+06,
0.174774686E+06, 0.175736558E+06, 0.176712352E+06, 0.177700040E+06, 0.178699106E+06,
0.179709512E+06, 0.180731620E+06, 0.181766307E+06, 0.182814343E+06, 0.183876525E+06,
0.184954817E+06, 0.186050290E+06, 0.187161233E+06, 0.188285994E+06, 0.189423534E+06,
0.190574069E+06, 0.191737635E+06, 0.192913719E+06, 0.194104807E+06, 0.195313945E+06,
0.196540986E+06, 0.197783113E+06, 0.199036424E+06, 0.200300179E+06, 0.201574293E+06,
0.202858665E+06, 0.204154096E+06, 0.205458732E+06, 0.206767566E+06, 0.208079174E+06,
0.209393060E+06, 0.210707027E+06, 0.212013409E+06, 0.213306648E+06, 0.214585746E+06,
0.215850982E+06, 0.217081712E+06, 0.218263024E+06, 0.219407918E+06, 0.220521581E+06,
0.221558365E+06, 0.222493116E+06, 0.223372462E+06, 0.224166053E+06, 0.224805399E+06,
0.225340329E+06, 0.225784503E+06, 0.225959322E+06, 0.226027519E+06, 0.225909648E+06,
0.225538722E+06, 0.225003136E+06, 0.224185757E+06, 0.223086575E+06, 0.221755830E+06,
0.220014998E+06, 0.218083899E+06, 0.215712952E+06, 0.213132640E+06, 0.210123945E+06,
0.206875823E+06, 0.203154809E+06, 0.199302180E+06])

um = 0.602393949E+00
up = 0.157520013E+00

um = -um / (2*pi)
up = -up / (2*pi)
dpsi = (um - up)/(MESHR-1)

my_pres = pres
my_fpol = fpol

for i in range(MESHR-1, 1):
  pprim_c = 0.5 * (pprime[i+1] + pprime[i])
  fprim_c = 0.5 * (ffprim[i+1] + ffprim[i])
  
  my_pres[i]=my_pres[i+1] + pprim_c*dpsi
  my_fpol[i]=sqrt(my_fpol[i+1]**2 + 2*fprim_c*dpsi)
  


psi = np.linspace(0, 1, len(fpol))

ax1= pyplot.subplot(231)
pyplot.scatter(psi, pres)
pyplot.ylabel('pres')
pyplot.xlabel('psi')
pyplot.grid(True)

ax = pyplot.subplot(232)
pyplot.scatter(psi, pprime)
pyplot.ylabel('pprim')
pyplot.xlabel('psi')
pyplot.grid(True)

ax = pyplot.subplot(233)
pyplot.scatter(psi, my_pres)
pyplot.ylabel('pres_der')
pyplot.xlabel('psi')
pyplot.grid(True)

ax = pyplot.subplot(234)
pyplot.scatter(psi, fpol)
pyplot.ylabel('fpol')
pyplot.xlabel('psi')
pyplot.grid(True)

ax = pyplot.subplot(235)
pyplot.scatter(psi, ffprim)
pyplot.ylabel('ffprim')
pyplot.xlabel('psi')
pyplot.grid(True)

ax = pyplot.subplot(236)
pyplot.scatter(psi, my_fpol)
pyplot.ylabel('fpol_der')
pyplot.xlabel('psi')
pyplot.grid(True)

pyplot.show()

# ax = pyplot.subplot(236)
# pyplot.scatter(psi[0:-1], ffprim[0:-1]/fpol_der)
# pyplot.ylabel('fpol multuplicator')
# pyplot.xlabel('psi')
# pyplot.grid(True)

# ax = pyplot.subplot(236)
# pyplot.scatter(psi[0:-1], pprim[0:-1]/pres_der)
# pyplot.ylabel('pres multuplicator')
# pyplot.xlabel('psi')
# pyplot.grid(True)

# pyplot.figure()
# pyplot.scatter(psi[0:-2], pprime[0: -2]/pres_der[1: :])
# pyplot.ylabel('pres multuplicator')
# pyplot.xlabel('psi')
# pyplot.grid(True)
# pyplot.savefig("Pics/pres_multuplicator_right.png", dpi=240)

# pyplot.figure()
# pyplot.scatter(psi[0:-1], (2*fpol[0:-1]*fpol_der)/ffprim[0:-1])
# pyplot.ylabel('fpol multuplicator')
# pyplot.xlabel('psi')
# pyplot.grid(True)
# pyplot.savefig("Pics/fpol_multuplicator.png", dpi=240)
# pyplot.show()

print(1)
# ffprim
# ffprim
# pres
# pprim
# from  write_eqdsk_equlx import *