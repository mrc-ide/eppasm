#############################################################################################################
## Working out how to get the detection rates going #########################################################
#############################################################################################################

det_I500A <- rep(0, 720)
det_I350A <- rep(0, 720) 
det_I200A <- rep(0, 720)
det_I0A <- rep(0, 720)

# annual rates, 1980-1985
det_I500A[181:431] <- 0
det_I350A[181:431] <- 0
det_I200A[181:252] <- 0.01
det_I0A[181:252] <- 0.05


### linear increase 1986 - end2000 ###
det_I200A[253:432] <- 0.01 + seq(1,180,1)*0.01
det_I0A[253:432] <- 0.05 + seq(1,180,1)*0.04


# # baseline rate - 1997 - end2000
det_I500A[385:432] <- 0.02 # baseline rate
det_I350A[385:432] <- 0.02 # baseline rate


### linear increase 2001-2009 ###
det_I500A[433:540] <- det_I500A[432] + seq(1,108,1)*0.02 
det_I350A[433:540] <- det_I350A[432] + seq(1,108,1)*0.02
det_I200A[433:540] <- det_I200A[432] + seq(1,108,1)*0.03
det_I0A[433:540] <- det_I0A[432] + seq(1,108,1)*0.05


### linear increase 2010-2016 ###
det_I500A[541:624] <- det_I500A[540] + seq(1,84,1)*0.021 
det_I350A[541:624] <- det_I350A[540] + seq(1,84,1)*0.021
det_I200A[541:624] <- det_I200A[540] + seq(1,84,1)*0.03
det_I0A[541:624] <- det_I0A[540] + seq(1,84,1)*0.05


### keep 2014 values for 2015-2025 ###
det_I500A[625:720] <- det_I500A[624]
det_I350A[625:720] <- det_I350A[624]
det_I200A[625:720] <- det_I200A[624]
det_I0A[625:720] <- det_I0A[624]

###################################################################################################
## So that works for single vectors of rates, lets try and scale it up to our arrays ##############
###################################################################################################

detec_array <- array(0, c(dim(brazil$fp$diagn_rate)))

## from 1980 change the rates for 3&4 and 5&6&7

detec_array[c(3:4),,,11:16] <- 0.01
detec_array[5:7,,,11:16] <- 0.05

### linear increase 1986 - 2000

detec_array[3:4,,,17:31] <- sweep(detec_array[3:4,,,17:31],4,(0.01 + seq(1,15,1) * 0.01),"+")
detec_array[5:7,,,17:31] <- sweep(detec_array[5:7,,,17:31],4,(0.05 + seq(1,15,1) * 0.04),"+")

detec_array[1,,,27:31] <- 0.02
detec_array[2,,,27:31] <- 0.02

## linear increase 2001 to 2009 ###

detec_array[1,,,32:40] <- sweep(detec_array[1,,,32:40], 3, (0.02 + seq(1,9,1) * 0.021), "+")
detec_array[2,,,32:40] <- sweep(detec_array[2,,,32:40], 3, (0.02 + seq(1,9,1) * 0.021), "+")
detec_array[3:4,,,32:40] <- sweep(detec_array[3:4,,,32:40], 4, (detec_array[3,1,1,31] + seq(1,9,1) * 0.03), "+")
detec_array[5:7,,,32:40] <- sweep(detec_array[5:7,,,32:40], 4, (detec_array[5,1,1,31] + seq(1,9,1) * 0.07), "+")

### linear increase 2010 to 2015 #####

detec_array[1,,,41:46] <- sweep(detec_array[1,,,41:46], 3, (detec_array[1,1,1,40] + seq(1,6,1) * 0.021), "+")
detec_array[2,,,41:46] <- sweep(detec_array[2,,,41:46], 3, (detec_array[2,1,1,40] + seq(1,6,1) * 0.021), "+")
detec_array[3:4,,,41:46] <- sweep(detec_array[3:4,,,41:46], 4,(detec_array[3,1,1,40] + seq(1,6,1) * 0.03), "+")
detec_array[5:7,,,41:46] <- sweep(detec_array[5:7,,,41:46], 4,(detec_array[5,1,1,40] + seq(1,6,1) * 0.11), "+")

### same values during proj period 

detec_array[1,,,47:52] <- detec_array[1,,,46] 
detec_array[2,,,47:52] <- detec_array[2,,,46]
detec_array[3:4,,,47:52] <- detec_array[3:4,,,46]
detec_array[5:7,,,47:52] <- detec_array[5:7,,,46]

matplot(detec_array[,,1,1],type = "l")
for(i in 2:52){
  matplot(detec_array[,,1,i], type = "l")
  Sys.sleep(1)
}


plot(detec_array[7,9,1,],type = "l",col="blue")
abline(v=46, col = "red")
abline(v=40,col="red")
abline(v=31, col="red")
abline(v=16, col="red")
abline(v= 11, col="red")
lines(detec_array[4,9,1,], type = "l", col="blueviolet")
lines(detec_array[2,9,1,], type = "l", col="forestgreen")
lines(detec_array[1,9,1,], type = "l", col="cyan")
