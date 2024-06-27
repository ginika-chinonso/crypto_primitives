# CIRCOM_GKR

Prove and verify any circom code with a GKR backend


## Steps
- Install the binary
- Write circom code
- Run binary with path to circom code



### Note:
This project is still under development and may contain bugs. It is not to be used in production until it has been audited.






// RUST_BACKTRACE=1 cargo run prove --input-path ./src/test_assets/test.circom --output-path ./hello.txt



// prove
// cargo run prove --circom-file-path test.circom --witness-path ./witness.json