# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2017-05-16 07:35
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('seqr', '0009_uploadedfileforfamily_uploadedfileforindividual'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='variantannotation',
            unique_together=set([]),
        ),
        migrations.AlterIndexTogether(
            name='variantannotation',
            index_together=set([]),
        ),
        migrations.RemoveField(
            model_name='variantannotation',
            name='created_by',
        ),
        migrations.RemoveField(
            model_name='variantnote',
            name='gene_id',
        ),
        migrations.RemoveField(
            model_name='variantnote',
            name='molecular_consequence',
        ),
        migrations.RemoveField(
            model_name='variantnote',
            name='transcript_id',
        ),
        migrations.AddField(
            model_name='variantnote',
            name='variant_annotation',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='variantnote',
            name='variant_genotypes',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='varianttag',
            name='variant_genotypes',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='varianttagtype',
            name='category',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.RemoveField(
            model_name='varianttag',
            name='variant_annotation',
        ),
        migrations.AddField(
            model_name='varianttag',
            name='variant_annotation',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.DeleteModel(
            name='VariantAnnotation',
        ),
    ]