# Music Unswinger
A program for adding or removing swing to music. Works by taking each beat of the song, splitting it in half and then stretching each half to fit swing time.

Eg. ``[0, 1/2, 1] -> [0, 2/3, 1]``

## Program usage:

```shell
./MusicUnswinger.exe [File] [Swing] [BPM] [Offset]
```

- **File**	- The song file you would like to be altered
- **Swing** - Must be either ``add`` or ``remove``
- **BPM**	- The BPM of the song
- **Offset** - The offset of the start of the song (in seconds). For example, 1.5s denotes the main downbeat of the song being on 

Although offset can be 0, it is still **mandatory**.

## Dependencies
- [libsndfile](https://github.com/libsndfile/libsndfile.git)
- [FTTW3](https://www.fftw.org/)