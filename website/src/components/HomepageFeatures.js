import React from 'react';
import clsx from 'clsx';
import styles from './HomepageFeatures.module.css';

const FeatureList = [
  {
    title: 'Supporting your “omics” research',
    Svg: require('../../static/img/red.svg').default,
    description: (
      <>
        WARP pipelines process a range of genomic, genotyping, transcriptomic and epigenomic data.
      </>
    ),
  },
  {
    title: 'Delivering production-level quality',
    Svg: require('../../static/img/yellow.svg').default,
    description: (
      <>
        Pipelines are used for production by the Broad Genomics Platform and the Human Cell Atlas; they are scientifically validated, scalable, and cost- and cloud- optimized
      </>
    ),
  },
  {
    title: 'Empowering reproducible science',
    Svg: require('../../static/img/blue.svg').default,
    description: (
      <>
        Pipelines are versioned and open-source, allowing you to track and share exactly how your data was processed.
      </>
    ),
  },
];

function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center">
        <Svg className={styles.featureSvg} alt={title} />
      </div>
      <div className="text--center padding-horiz--md">
        <h3>{title}</h3>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
