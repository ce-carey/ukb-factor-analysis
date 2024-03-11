set -ex

cd ../../lavaan_robust_scaled
tar -zcvf lavaan.tar.gz lavaan
cd ../TRAINING/redux10point2

zone=us-east1-b
prefix=ccarey
script_name=cfa

name=$prefix-$script_name-training
gcloud compute instances create --zone $zone --machine-type m1-ultramem-80 --boot-disk-size=1500GB --image cfa-benchmarking $name &

wait

sleep 5 # give the instances a moment to fully come alive

gcloud compute scp --zone $zone $script_name.r ../../lavaan_robust_scaled/lavaan.tar.gz ../../imp1.Rds ../train_withdraw_outliers.txt ../../full.txt muthen1984_fullweight.Rds lav_model_estimate.Rds $name:. 
