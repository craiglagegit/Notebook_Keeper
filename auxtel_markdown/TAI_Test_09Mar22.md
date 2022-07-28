---
jupyter:
  jupytext:
    formats: ipynb,markdown//md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.0
  kernelspec:
    display_name: LSST
    language: python
    name: lsst
---

## TAI testing - NCSA
Craig Lage - Mar 9, 2022 \
Repeating - May 4, 2022

```python
from astropy.time import Time
t = Time('2022-03-08', scale='utc')
t.unix_tai - t.unix
```

```python
from astropy import __version__ as astropyVersion
print(astropyVersion)
```

```python

```
