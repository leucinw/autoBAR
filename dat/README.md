
# README

| Filename | Remark |
|---|---|
|gas.key | This exact file is not used in jobs unless no key provided in the respective `settings.yaml`. |
|liquid.key | This exact file is not used in jobs unless no key provided in the respective `settings.yaml`. |
|orderparams_courser | This exact file is used in jobs if courser lambda windows enabled. If user provided another lambda_window_file in `settings.yaml`, the user file will be used instead.|
|orderparams_default | This exact file is used in jobs if default lambda windows enabled. If user provided another lambda_window_file in `settings.yaml`, the user file will be used instead.|
|settings.yaml | This exact file is just an example. The file in individual `targets` directory is the one actually used. |
|tinker.env | This exact file is always used in jobs. |



