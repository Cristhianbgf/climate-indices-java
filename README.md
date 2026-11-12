# climind

Drought and climate extreme indices in Java, checked against [xclim](https://github.com/Ouranosinc/xclim). The SPI comes first; SPEI and the ETCCDI indices are planned.

I could not find a maintained Java library for these indices. People who work in Java (THREDDS, ERDDAP, GeoServer or institutional systems) usually call Python or R to get them. The aim here is to compute them in Java and to document how the results compare with the reference implementations.

Work in progress. See the open issues for what comes next.

## Build

Requires Java 21 and Maven.

```bash
mvn test
```

## Licence

MIT. See [LICENSE](LICENSE).
