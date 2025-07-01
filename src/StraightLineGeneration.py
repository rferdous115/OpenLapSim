# Generating a straight‐line track file in the exact two‐column format
distances = range(0, 6951, 10)
with open('trackFiles/straight_track.txt', 'w') as f:
    for d in distances:
        f.write(f"{d}\t0\n")

# Show a preview
with open('trackFiles/straight_track.txt') as f:
    print(''.join([next(f) for _ in range(10)]))