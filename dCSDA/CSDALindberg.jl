#=
Stand alone program to obtain CSD using the Lindberg Opeartors
This program works with files that have been already pre procesed, that is,
the failling electrodes have been singled out and the physiologically active
channels have also been tagged.
=#



using HDF5
push!(LOAD_PATH, ".")
# Auxiliary files with operators and conversion tools
using LindbergOperators, ParaSets


if length(ARGS)==0
    error("A hdf5 data file is needed to work.")
else
    archivo=ARGS[1]
    println("Working with the file: ", archivo)
end

arx=h5open(archivo)

#Load LFP data and lists of failing and active electrodes.

# Please check the name that you used for saving the LFP

LFPName="LFPSaturados" #example
LFP=read(arx[LFPName])

saturadosarray=read(arx["CanalesSaturados"])
respuestasarray=read(arx["Canalesrespuesta"]);
frecuencia=read(arx["freq"])
close(arx)

saturados=arraytoset(saturadosarray)
respuestas=arraytoset(respuestasarray)


# Copy the data on an auxiliary Array
lfpParchado=copy(LFP)

yes

#Set to zero the unusable channels
for m in saturados
    q=m[1]
    p=m[2]
    lfpParchado[q,p,:]=0
end


#
# Edges are irrelevant
listaredux=TiraOrillas(saturados)


# Smooth artificial zeros
for m in listaredux
        q=m[1]
        p=m[2]
        vecinos=vecindad8(m)
        lfpParchado[q,p,:]=promediasobreconjunto(vecinos,lfpParchado)
end



(mu,nu,lu)=size(lfpParchado)


# Gaussian Temporal Smoothing
lfpplanchado=zeros(mu,nu,lu)
for j=1:mu,l=1:nu
    porromponpon=vec(lfpParchado[j,l,:])
    lfpplanchado[j,l,:]=GaussSuavizarTemporal(porromponpon)
end


aux1=zeros(mu,nu,lu)
aux2=zeros(mu,nu,lu)
# Gaussian Smoothing and dCSDA
#Posteriormente sacamos el dCSD.
for t=1:lu
    aux1[:,:,t]=GaussianSmooth(lfpplanchado[:,:,t])
    aux2[:,:,t]=DiscreteLaplacian(aux1[:,:,t])
end
CSD=-aux2


# we save in the same file the  new data.
h5open(archivo, "r+") do file
    write(file, "CSDALindberg", CSD)  # alternatively, say "@write file A"
end


