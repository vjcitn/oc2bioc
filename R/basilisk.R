
ocbenv <- basilisk::BasiliskEnvironment(envname="ocbenv",
    pkgname="oc2bioc", packages=c("python=3.8.2", "requests==2.22.0", 
      "pytz==2019.3"),
    pip=c("open-cravat==2.2.1"))
