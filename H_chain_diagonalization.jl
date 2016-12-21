
using PyPlot


## ========== number of particles basis ================


n = 10
function mysub2ind(inds::Vector{Int64},n::Int64)
  res = inds[n]-1
  for i=n-1:-1:1
    res *= 2
    res += inds[i]-1
  end
  res+1
end

si1 = (fill(2,n)...)
#inds = [ind2sub(si1,5)...]
#mysub2ind(inds,10)


row = Vector{Int64}(0)
for i = 1:2^n
  inds = [ind2sub(si1,i)...]
  value = sum(inds)
  #ind = Array(inds)
  if value == (3*n/2)
    push!(row,i)
  end
end


Hrow = Vector{Int64}(0)
Hcol = Vector{Int64}(0)
Hval = Vector{Float64}(0)
for i = 1:size(row)[1]
  row_inds = [ind2sub(si1,row[i])...]
  for j = 1:size(row)[1]
    col_inds = [ind2sub(si1,row[j])...]
      diff = row_inds - col_inds
      diffnz = find(diff)
      #if size(diffnz)[1] == 2 && ((diffnz[2]-diffnz[1] == 1) || (diffnz[2]-diffnz[1] == 2))
      if size(diffnz)[1] == 2
        @show row_inds
        if (diffnz[2]-diffnz[1] == 1)
            #@show diffnz,size(diffnz)
            @show 1,row_inds,col_inds,i,j
            v = 1
            push!(Hrow,i)
            push!(Hcol,j)
            push!(Hval,v)
        end
        if (diffnz[2]-diffnz[1] == 2)
            #@show diffnz,size(diffnz)
            @show 2,row_inds,col_inds,i,j,diffnz[1]
            if row_inds[diffnz[1]+1] == 1
              v = 1
            else
              v = -1
            end
            push!(Hrow,i)
            push!(Hcol,j)
            push!(Hval,v)
        end
      end
  end
end
H = sparse(Hrow,Hcol,Hval)
flush(STDOUT)
Hamiltonian = full(H)
f1 = open("/Users/yifantian/Desktop/Tensor/Hubbard/H chain/Hamiltonian2.txt","w")
for i = 1:size(row)[1]
  for j = 1:size(row)[1]
    s = @sprintf "%2.1d " H[i,j]
    write(f1,s)
  end
  write(f1,"\n")
end
close(f1)

#=
Hamiltonian == Hamiltonian1 + Hamiltonian2
evn = eigs(Hamiltonian1;nev=1,which=:SR)
energy = evn[1]
evn = eigs(Hamiltonian2;nev=1,which=:SR)
energy = evn[1]
evn = eigs(Hamiltonian;nev=1,which=:SR)
energy = evn[1]
=#
evn = eigs(H;nev=1,which=:SR)
energy = evn[1]

@show evn[2]


Cmn_1 = zeros(9)
Cmn_2 = zeros(8)
#m = 2
for m = 1:9
#for n = 3:9
  #Cmn = 0.0
  for i = 1:size(evn[2])[1]
    for j = 1:size(evn[2])[1]
      row_inds = [ind2sub(si1,row[i])...]
      col_inds = [ind2sub(si1,row[j])...]
      diff = row_inds - col_inds
      diffnz = find(diff)
      if (size(diffnz)[1] == 2)
        #if (diff[m] == 1 && diff[n] == -1)
        if (diff[m] == 1 && diff[m+1] == -1)
          #@show diff,i,j,evn[2][i],evn[2][j]
          Cmn_1[m] += (-1*evn[2][i]*evn[2][j])
          #Cmn += (1*evn[2][i]*evn[2][j])
        end
          if (m < 9)
            if (diff[m] == 1 && diff[m+2] == -1)
              #@show diff,i,j,evn[2][i],evn[2][j]
              #Cmn_2[m] += (1*evn[2][i]*evn[2][j])
              #Cmn += (1*evn[2][i]*evn[2][j])
              if row_inds[diffnz[1]+1] == 1
                Cmn_2[m] += (-1*evn[2][i]*evn[2][j])
              else
                Cmn_2[m] += (1*evn[2][i]*evn[2][j])
              end
            end
          end
      end
    end
  end
#@show Cmn
end
#@show Cmn
E_difference = energy - (-2*sum(Cmn_1)) - (-2*sum(Cmn_2))


###＝＝＝＝＝＝＝ new way concerning anticommutation

## ====== Cdag_i C_i+2 ====
Cmn = 0.0
m = 1                                             # (1,3) term
for i = 1:size(evn[2])[1]
  for j = 1:size(evn[2])[1]
    row_inds = [ind2sub(si1,row[i])...]
    col_inds = [ind2sub(si1,row[j])...]
    diff = row_inds - col_inds
    diffnz = find(diff)
    if (size(diffnz)[1] == 2)                     # hopping term
      #if (diff[m] == 1 && diff[n] == -1)
      if (diff[m] == 1 && diff[m+2] == -1)
        @show row_inds,col_inds,diff,i,j,evn[2][i],evn[2][j]
        Cmn += evn[2][i]*evn[2][j]*(1)
        #Cmn += (1*evn[2][i]*evn[2][j])
      end
    end
  end
end
Cmn
## ====== Cdag_i C_i+3 ====







## ======

Cmn = zeros(9)
m = 1
#for m = 1:9
for n = 2:10
  #Cmn = 0.0
  for i = 1:size(evn[2])[1]
    for j = 1:size(evn[2])[1]
      row_inds = [ind2sub(si1,row[i])...]
      col_inds = [ind2sub(si1,row[j])...]
      diff = row_inds - col_inds
      diffnz = find(diff)
      if (size(diffnz)[1] == 2)
        if (diff[m] == 1 && diff[n] == -1)
        #if (diff[m] == 1 && diff[m+1] == -1)
          #@show diff,i,j,evn[2][i],evn[2][j]
          #Cmn_1[m] += (1*evn[2][i]*evn[2][j])
          Cmn[n-1] += (1*evn[2][i]*evn[2][j])
        end
      end
    end
  end
#@show Cmn
end

x = [2:1:10]
plot(x, Cmn, color="red", linewidth=2.0, linestyle="--")
title("<Cdagup_i*Cup_j>, N=10 lattice")

x = [1:1:9]
plot(x, Cmn_1, color="red", linewidth=2.0, linestyle="--")
title("<Cdagup_i*Cup_i+1>, N=10 lattice")

#####
Nsites = zeros(10)
m = 5
for i = 1:size(evn[2])[1]
  for j = 1:size(evn[2])[1]
    row_inds = [ind2sub(si1,row[i])...]
    col_inds = [ind2sub(si1,row[j])...]
    diff = row_inds - col_inds
    diffnz = find(diff)
    if (size(diffnz)[1] == 0)
      for l = 1:10
        Nsites[l] += evn[2][i]*evn[2][j]*(row_inds[l]-1)
      end
    end
  end
end
@show Nsites



##======== energy levels basis ==============

t = 1
U = 0
n = 10

row = Vector{Int64}(0)
col = Vector{Int64}(0)
val = Vector{Float64}(0)

function dopush!(i,j,v)      #dynamically allocate
  push!(row,i)
  push!(col,j)
  push!(val,v)
end

for i=1:n-1
  #dopush!(i,i,U)
  jvalue1 = (i+1)
  jvalue2 = (i+1)
  #if (i == 7)
  #  jvalue1 = 8
  #elseif (i == 1)
  #  jvalue2 = 8
  #end
  dopush!(i,jvalue1,-t)
  dopush!(jvalue2,i,-t)
end

for i=1:n-2
  #dopush!(i,i,U)
  jvalue1 = (i+2)
  jvalue2 = (i+2)
  #if (i == 7)
  #  jvalue1 = 8
  #elseif (i == 1)
  #  jvalue2 = 8
  #end
  dopush!(i,jvalue1,-t)
  dopush!(jvalue2,i,-t)
end

H = sparse(row,col,val)
flush(STDOUT)

full(H)

#evn = eigs(H;nev=1,which=:SR)
#evn[1]

evn = eig(full(H))
Elevels = evn[1]
Energy = sum(Elevels[1:5])

#evn = eigs(H,nev = 10)
#Elevels = sort(evn[1])
#sum(Elevels[1:5])

H12 = zeros(10,10)
H12[1,2] = 1

H36 = zeros(10,10)
H36[3,6] = 1

H = zeros(10,10)
H13[1,3] = 1

Corr = cell(6)
for n = 5:10
  Cij = 0.0
  cMatrix = zeros(10,10)
  cMatrix[5,n] = 1
  for i = 1:5
    Cij += transpose(evn[2][:,i])*cMatrix*evn[2][:,i]
  end
  Corr[n-4] = Cij[1]
end

x = [5:1:10]
plot(x, Corr, color="red", linewidth=2.0, linestyle="--")
title("<Cdagup_i*Cup_j>, N=10 lattice")


C12 = 0.0
for i = 1:5
  C12 += transpose(evn[2][:,i])*H12*evn[2][:,i]
end
@show C12




















####
