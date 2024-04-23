# SPACE EFFICIENT BINARY MERKLE TREE IMPLEMENTATION

This file contains an implementation of a space efficient binary merkle tree

NOTE: 

Normally, it is expected that the length of the leaves of a binary merkle tree is a power of 2. This implementation is not space efficient as it pads the leaves with dummy values up to the next power of 2.

This implementation doesnt require that restriction but only requires that the number of leaves are even thereby cutting down on the space consumption.