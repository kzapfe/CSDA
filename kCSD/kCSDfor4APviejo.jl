using JLD

nombre="../Actividad4AP/DatosLFP4APusadosPaper.jld"
archivo=load(nombre)
LFP=archivo["LFP"];

K=readdlm("K_Dura_r_0.25-4095.dat");
Ktilde=readdlm("KTilde_Dura_r_0.25-4095.dat");


Kinverso=inv(K);
KTTKInv=transpose(Ktilde)*Kinverso;

tmax=size(LFP, 3)
CSD=zeros(LFP)

#t=200

for t=1:tmax
    testLFP=LFP[:,:,t]
    testLFPvector1=reshape(transpose(testLFP),4096);
    tlfpred=testLFPvector1[2:end];
    CSDtentativa=zeros(4096);
    CSDtentativa[2:end]=transpose(KTTKInv)*tlfpred;
    CSD[:,:,t]=reshape(CSDtentativa,64,64)
end





paguardar=load(nombre)
paguardar["kCSDCorrecta"]=CSD
save(nombre,paguardar)
